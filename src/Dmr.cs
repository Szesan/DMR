// The full DMR model: reconstruction, exon validation, and swap-and-amplify. See README §7.3-7.5.
using System;
using System.Collections.Generic;
using System.Linq;
/// <summary>The generic DNA mechanism of the DMR model: the framework's processes made executable.</summary>
public class Dmr : IPreservationModel
{
    readonly IGenome genome;

    public IGenome Genome => genome;

    public Dmr(IGenome genome)
    {
        this.genome = genome;
    }

    /// <summary>Validates repeating genomic units within a tandem array</summary>
    public Pool<ISequence> ConsensusOnTandemArray(DoubleStrand segmentsDelimitedBy5hmuMotifs)
    {
        (ISequence senseStrand, ISequence template) = Enzyme.OpenDoubleHelix(segmentsDelimitedBy5hmuMotifs);

        Pool<ISequence> validatedSegments = new Pool<ISequence>();

        template.NewRead();

        while (template.HasNextBase())
        {
            template.ReadNextBase();   // discard the segment's leading 5hmu border
            
            RNA segmentProbe = Enzyme.TranscribeUntilBorder(template, genome);

            segmentProbe.Prepend(Base.Border);  

            template.Cleave();
            if (Enzyme.FindsSeatingOn(segmentProbe, template) && template.Base == Base.Border)
            {
                Enzyme.ReleaseSequence(segmentProbe);
                validatedSegments.Add(segmentProbe);
            }

            template.NewRead();                                   
        }

        return validatedSegments;
    }

    /// <summary>Validates the exons that tile up to a certified gene</summary>
    public Pool<ISequence> ValidateExonsFromGene(ISequence gene)
    {
        ISequence template = Enzyme.SynthesizeComplementaryStrand(gene, genome).ComplementaryStrand;

        Pool<ISequence> rootPool = Genome.FetchSequencesByFlankMotifs(Settings.ModulonLeftFlank, Settings.ModulonRightFlank);

        Pool<ISequence> bodies = new Pool<ISequence>();

        while (rootPool.IsNotEmpty())
        {
            ISequence body = rootPool.PullOne();

            body.NewRead();
            while (body.Base != Base.Border) body.ReadNextBase();   // walk off the [modulon] tag
            body.ReadNextBase();                                    // step past the internal border
            body.Cleave();                                          // discard [modulon] h -> pure exon body

            bodies.Add(body);
        }

        ISequence senseStrand = template.CreateEmptyBackbone();
        Pool<ISequence> validated = new Pool<ISequence>();

        while (!senseStrand.ContainsNoneOf(Base.Empty) && bodies.IsNotEmpty())
        {
            ISequence exon = bodies.PullOne();

            bool match;
            while ((match = Enzyme.FindsSeatingOn(exon, template)) && !exon.ClosesCleanlyAgainst(senseStrand))
                Enzyme.Shift(exon);

            if (match)
            {
                Enzyme.Bond(senseStrand, exon);
                Enzyme.ReleaseSequence(exon);
                validated.Add(exon);
            }
        }

        return validated;
    }

    /// <summary>Lays the segments out head-to-tail as one strand, each capped by a 5hmu boundary mark</summary>
    public DoubleStrand LayDownTandemArray(Pool<ISequence> segments)
    {
        ISequence tandemArray = Genome.CreateNewSequence();

        tandemArray.AddBase(Base.Border);     

        while(segments.IsNotEmpty())
        {
            ISequence segment = segments.PullOne();

            segment.NewRead();
            
            if (!segment.HasNextBase()) continue;

            while (segment.HasNextBase())
            {
                tandemArray.AddBase(segment.Base);
                segment.ReadNextBase();
            }

            tandemArray.AddBase(Base.Border);
        }            

        return Enzyme.SynthesizeComplementaryStrand(tandemArray, genome);
    }

    /// <summary>Builds [Word][Modulon][Boundary] decoder DNA</summary>
    public Pool<DoubleStrand> BuildRecipeDecoders(Pool<ISequence> maps, Pool<ISequence> modulonTaggedRoots)
    {
        Pool<DoubleStrand> decoders = new Pool<DoubleStrand>();
        Pool<ISequence> modulonPool = new Pool<ISequence>();

        while(modulonTaggedRoots.IsNotEmpty())
        {
            ISequence root = modulonTaggedRoots.PullOne();
            root.NewRead();

            ISequence modulonTemplate = Enzyme.TranscribeUntilBorder(root, genome).Sequence; 
            modulonPool.Add(modulonTemplate);
        }

        Pool<ISequence> enrichedModulons = Enzyme.Enrich(modulonPool, Settings.EnrichmentAbundant, genome);

        while (maps.IsNotEmpty())
        {
            ISequence map = maps.PullOne();

            while (enrichedModulons.IsNotEmpty())
            {
                ISequence modulon = enrichedModulons.PullOne();

                if (!Enzyme.FindsSeatingOn(modulon, map)) continue;

                modulon.AddBase(Base.Border);

                map.AddBase(Base.Border);

                decoders.Add(new DoubleStrand(map, modulon));

                break; 
            }
        }

        return decoders;
    }

    public DoubleStrand DecodeRecipeToModulonSequence(DoubleStrand recipe, Pool<DoubleStrand> decoders)
    {
        DoubleStrand recipeTemplate = new DoubleStrand(recipe.SenseStrand.CreateEmptyBackbone(), recipe.ComplementaryStrand);

        Pool<DoubleStrand> pool = Enzyme.Enrich(decoders, Settings.EnrichmentAbundant, genome);

        while(pool.IsNotEmpty())
        {
            DoubleStrand decoder = pool.PullOne();

            recipeTemplate.NewRead();

            while(recipeTemplate.HasNextBase())
            {
                //Actually double-doubleStrand, treating as such simplifies lockstep
                DoubleStrand decodedRecipe = new DoubleStrand(decoder,recipeTemplate);

                decodedRecipe.NewRead();

                while(decoder.ClosesAtHead(recipeTemplate)
                    && decoder.ComplementaryStrand.Base == Base.Empty) decodedRecipe.ReadNextBase();

                if(decoder.ComplementaryStrand.Base != Base.Empty)
                {
                    decodedRecipe.NewRead();

                    while(decoder.SenseStrand.Base == Base.Empty) decodedRecipe.ReadNextBase();

                    DoubleStrand leftFlank = (DoubleStrand)(recipeTemplate.Cleave());   // guide-so-far, kept
                    decoder.Cleave();                                          // drop probe's leading e's

                    decodedRecipe.NewRead();

                    while(decoder.ComplementaryStrand.Base == Base.Empty) decodedRecipe.ReadNextBase();

                    decodedRecipe.Cleave();                                          // drop [W], keep [M][S]

                    recipeTemplate = (DoubleStrand)(leftFlank.Ligate(decoder).Ligate(recipeTemplate));

                    if(recipeTemplate.SenseStrand.ContainsNoneOf(Base.Empty)) return recipeTemplate; 

                    break;
                }
                else Enzyme.Shift(decoder);
            }
        } 
        
        return recipeTemplate;    
    }

    public Pool<ISequence> AssembleGeneVariants(ISequence guideSense, Pool<ISequence> rootPool)
    {
        Pool<ISequence> variants = new Pool<ISequence>();
        Pool<ISequence> rejected = new Pool<ISequence>();
        ISequence template = Enzyme.SynthesizeComplementaryStrand(guideSense, genome).ComplementaryStrand;

        while (rootPool.IsNotEmpty())
        {
            ISequence root      = rootPool.PullOne();
            ISequence rootProbe = Enzyme.CopyStrand(root, genome);

            template.NewRead();
            while (template.HasNextBase())
            {
                DoubleStrand gene = new DoubleStrand(rootProbe, template);
                gene.NewRead();

                while (template.HasNextBase() && Enzyme.ClosesAtHead(rootProbe, template) && rootProbe.Base != Base.Border) gene.ReadNextBase();

                if (rootProbe.Base == Base.Border && template.Base == Base.Border)
                {
                    gene.NewRead();
                    while (rootProbe.Base == Base.Empty) gene.ReadNextBase();

                    ISequence leftFlank = (gene.Cleave() as DoubleStrand).ComplementaryStrand;

                    while (rootProbe.Base != Base.Border) gene.ReadNextBase();   // walk the probe's [M] to the Border
                    gene.ReadNextBase();                                     // skip the Border ([S] spacer)
                    gene.Cleave();                                           // probe remainder is the body [R]

                    template = leftFlank.Ligate(rootProbe).Ligate(template);
                    break;
                }
                else Enzyme.Shift(rootProbe);
            }

            if (!template.HasNextBase()) rejected.Add(root);                 // probe seated nowhere

            if (template.ContainsNoneOf(Base.Border) && template.ContainsNoneOf(Base.Empty))   // every slot a body -> variant complete
            {
                variants.Add(template);
                rootPool.Recycle(rejected);   
                template = Enzyme.SynthesizeComplementaryStrand(guideSense, genome).ComplementaryStrand;
            }
        }

        return variants;
    }

    public virtual List<ISequence> RunCycle(int time)
    {
        Pool<ISequence> genes = new Pool<ISequence>();
        Pool<ISequence> maps = new Pool<ISequence>();
        Pool<ISequence> recipes = new Pool<ISequence>();
 
        if (time.DueForMaintenance() || time.DueForReport())
        {
            maps = Genome.FetchSequencesByFlankMotifs(Settings.MapLeftFlank, Settings.MapRightFlank);
            maps = ConsensusOnTandemArray(LayDownTandemArray(maps));
            maps = ConsensusOnTandemArray(LayDownTandemArray(maps));
 
            recipes = Genome.FetchSequencesByFlankMotifs(Settings.RecipeLeftFlank, Settings.RecipeRightFlank);
            recipes = ConsensusOnTandemArray(LayDownTandemArray(recipes));
            recipes = ConsensusOnTandemArray(LayDownTandemArray(recipes));
 
            Pool<ISequence> decoderRoots = Genome.FetchSequencesByFlankMotifs(Settings.ModulonLeftFlank, Settings.ModulonRightFlank);
            Pool<DoubleStrand> decoders = BuildRecipeDecoders(new Pool<ISequence>(maps), decoderRoots);   // copy: maps must survive the refresh
 
            Func<ISequence, Pool<ISequence>> ReconstructGeneFromRecipe = recipe =>
            {
                DoubleStrand recipeDuplex = Enzyme.SynthesizeComplementaryStrand(recipe, genome);
                DoubleStrand guide = DecodeRecipeToModulonSequence(recipeDuplex, decoders);
 
                // Idiosyncratic damage: try reshuffled seatings until one ordering's consensus certifies (shared health).
                Pool<ISequence> rootPool = Genome.FetchSequencesByFlankMotifs(Settings.ModulonLeftFlank, Settings.ModulonRightFlank);
                Pool<Pool<ISequence>> orderings = rootPool.ShuffledOrderings(Settings.EnrichmentAbundant*Settings.EnrichmentAbundant, new Random(Settings.ShuffleSeed));
 
                Pool<ISequence> certified = new Pool<ISequence>();

                while (certified.IsEmpty() && orderings.IsNotEmpty())
                {
                    Pool<ISequence> variants  = AssembleGeneVariants(guide.SenseStrand, orderings.PullOne());
                    Pool<ISequence> singlyValidated = ConsensusOnTandemArray(LayDownTandemArray(variants));
                    certified = ConsensusOnTandemArray(LayDownTandemArray(singlyValidated));
                }
                return certified;
            };
 
            Pool<Pool<ISequence>> certifiedGenes = Genome.ProcessAtScale(ReconstructGeneFromRecipe, recipes);
 
            while (certifiedGenes.IsNotEmpty()) genes.Recycle(certifiedGenes.PullOne());
        }
 
        if (time.DueForMaintenance())
        {
            SwapAndAmplify(genes, maps, recipes);
        }
 
        return genes.ToList();   
    }

    public virtual void SwapAndAmplify(Pool<ISequence> genes, Pool<ISequence> maps, Pool<ISequence> recipes)
    {
        Pool<Pool<ISequence>> perGene = Genome.ProcessAtScale(ValidateExonsFromGene, genes);
        Pool<ISequence> validatedDuplexes = new Pool<ISequence>();
        while (perGene.IsNotEmpty())
        {
            Pool<ISequence> bodies = perGene.PullOne();
            while (bodies.IsNotEmpty())
                validatedDuplexes.Add(Enzyme.SynthesizeComplementaryStrand(bodies.PullOne(), genome));
        }

        Pool<ISequence> nextGenerationRoots = new Pool<ISequence>();
        Pool<ISequence> retiredModulons = new Pool<ISequence>();   // pre-advance tags: retire these lineages

        Pool<ISequence> roots   = Genome.FetchSequencesByFlankMotifs(Settings.ModulonLeftFlank, Settings.ModulonRightFlank);

        while (roots.IsNotEmpty())
        {
            ISequence exonBody = Enzyme.CopyStrand(roots.PullOne(), genome);

            exonBody.NewRead();
            while (exonBody.Base != Base.Border) exonBody.ReadNextBase();   // walk off the [modulon] tag
            ISequence modulonTag = exonBody.Cleave();                       // the [modulon] tag, pre-advance
            exonBody.ReadNextBase();                                        // step past the internal border
            exonBody.Cleave();                                              // now the pure exon body (sense)

            DoubleStrand bodyDuplex = Enzyme.SynthesizeComplementaryStrand(exonBody, genome);

            bool matches = false;
            Pool<ISequence> candidates = new Pool<ISequence>(validatedDuplexes);   // test against all, without draining the set
            while (candidates.IsNotEmpty() && !matches)
            {
                ISequence validated = candidates.PullOne();
                validated.NewRead();
                bodyDuplex.NewRead();
                matches = bodyDuplex.ClosesCleanlyAgainst(validated);
            }

            if (!matches) continue;   // exon body doesn't validate -> not a candidate root

            // exon body alone can't catch a mutated tag: the tag itself must appear in a certified map too
            ISequence modulonComplement = Enzyme.SynthesizeComplementaryStrand(modulonTag, genome).ComplementaryStrand;

            bool modulonMatches = false;
            Pool<ISequence> mapCandidates = new Pool<ISequence>(maps);
            while (mapCandidates.IsNotEmpty() && !modulonMatches)
            {
                ISequence map = mapCandidates.PullOne();
                ISequence modulonProbe = Enzyme.CopyStrand(modulonComplement, genome);
                map.NewRead();
                modulonProbe.NewRead();
                modulonMatches = Enzyme.FindsSeatingOn(modulonProbe, map);
            }

            if (!modulonMatches) continue;   // tag mutated beyond recognition -> pass on this root

            retiredModulons.Add(Enzyme.CopyStrand(modulonTag, genome));   // the pre-advance tag retires this lineage

            ISequence advancedModulon = modulonTag;
            Enzyme.AdvanceUnderStencil(advancedModulon, Enzyme.Stencil(Settings.ModulonStencil, genome));

            ISequence border = genome.CreateNewSequence();
            border.AddBase(Base.Border);
            nextGenerationRoots.Add(advancedModulon.Ligate(border).Ligate(exonBody));   // [advanced-modulon] h [exon]

        }

        // Advance each validated map's modulon under the map stencil; the archive word is shouldered.
        Pool<ISequence> newMaps = new Pool<ISequence>();

        while (maps.IsNotEmpty())
        {
            ISequence map = maps.PullOne();
            Enzyme.AdvanceUnderStencil(map, Enzyme.Stencil(Settings.MapStencil, genome));
            newMaps.Add(map);
        }

        Genome.EliminateByMotif(Settings.MapLeftFlank,     Settings.MapRightFlank);
        Genome.EliminateByMotif(Settings.RecipeLeftFlank,  Settings.RecipeRightFlank);

        while (retiredModulons.IsNotEmpty())
            Genome.EliminateByMotif(Settings.ModulonLeftFlank + retiredModulons.PullOne().ToString(), Settings.ModulonRightFlank);

        Pool<ISequence> newGenRoots = Enzyme.Balance(nextGenerationRoots, Settings.RootEnrichmentRate, genome);
        Pool<ISequence> newGenMaps = Enzyme.Balance(newMaps, Settings.IdentifierEnrichmentRate, genome);
        Pool<ISequence> newGenRecipes = Enzyme.Balance(recipes,   Settings.IdentifierEnrichmentRate, genome);

        Genome.NewRead();
        Genome.StoreDistributed(newGenRoots, Settings.ModulonLeftFlank, Settings.ModulonRightFlank);
        Genome.StoreDistributed(newGenMaps, Settings.MapLeftFlank,     Settings.MapRightFlank);
        Genome.StoreDistributed(newGenRecipes, Settings.RecipeLeftFlank,  Settings.RecipeRightFlank);
    }

}

