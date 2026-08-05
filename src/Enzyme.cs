// The molecular operations available to a model. See README §2.
using System;
using System.Linq;
public static class Enzyme
{
    /// <summary>Unzips the duplex into its two strands.</summary>
    public static (ISequence senseStrand, ISequence complementaryStrand) OpenDoubleHelix(DoubleStrand sequence)
    {
        return (sequence.SenseStrand, sequence.ComplementaryStrand);
    }

    /// <summary>Performs DNA replication/synthesis by assembling its complementary strand using Watson-Crick base pairing.</summary>
    public static DoubleStrand SynthesizeComplementaryStrand(ISequence templateStrand, IGenome genome)
    {
        ISequence complementaryStrand = genome.CreateNewSequence();

        templateStrand.NewRead();
        while (templateStrand.HasNextBase())
        {
            complementaryStrand.AddBase(GetWatsonAndCrickPair(templateStrand.Base));
            templateStrand.ReadNextBase();
        }

        return new DoubleStrand(templateStrand, complementaryStrand);
    }

    /// <summary>Replicates a strand verbatim by Watson-Crick synthesis: the complement of its complement.</summary>
    public static ISequence CopyStrand(ISequence strand, IGenome genome)
    {
        ISequence complement = SynthesizeComplementaryStrand(strand, genome).ComplementaryStrand;
        return SynthesizeComplementaryStrand(complement, genome).ComplementaryStrand;
    }

    /// <summary>Tests whether two opposing bases hybridize. An empty AP site pairs permissively.</summary>
    public static bool BasePairsWith(Base senseBase, Base complementaryBase)
    {
        if (senseBase == Base.Empty) return true;
        if (complementaryBase == Base.Empty) return true;

        if (GetWatsonAndCrickPair(senseBase) == complementaryBase) return true;

        return false;
    }

    /// <summary>Single-strand closure test AT THE CURRENT HEADS only.</summary>
    public static bool ClosesAtHead(ISequence probe, ISequence template)
    {
        if (!probe.HasNextBase() || !template.HasNextBase()) return false;
        return BasePairsWith(probe.Base, template.Base);
    }

    /// <summary>Returns the Watson-Crick complement of a base. (Empty - Anything, Border - Border)</summary>
    public static Base GetWatsonAndCrickPair(Base baseToPair)
    {
        if (baseToPair == Base.A) return Base.T;
        if (baseToPair == Base.T) return Base.A;
        if (baseToPair == Base.C) return Base.G;
        if (baseToPair == Base.G) return Base.C;
        if (baseToPair == Base.Border) return Base.Border;
        return Base.Empty;
    }

    /// <summary>Covalently bonds bases from a donor sequence into the empty AP sites of a host sugar-phosphate backbone.</summary>
    public static ISequence Bond(ISequence hostStrand, ISequence donorStrand)
    {      
        hostStrand.NewRead();
        donorStrand.NewRead();

        while(donorStrand.HasNextBase() && hostStrand.HasNextBase())
        {
            if(hostStrand.Base == Base.Empty) hostStrand.Base = donorStrand.Base;
            hostStrand.ReadNextBase();
            donorStrand.ReadNextBase();
        }

        return hostStrand;
    }

    public static DoubleStrand Bond(DoubleStrand hostStrand, DoubleStrand donorStrand)
    {
        hostStrand.NewRead();
        donorStrand.NewRead();

        while(donorStrand.HasNextBase() && hostStrand.HasNextBase())
        {
            if(hostStrand.SenseStrand.Base == Base.Empty) hostStrand.SenseStrand.Base = donorStrand.SenseStrand.Base;
            if(hostStrand.ComplementaryStrand.Base == Base.Empty) hostStrand.ComplementaryStrand.Base = donorStrand.ComplementaryStrand.Base;
            hostStrand.ReadNextBase();
            donorStrand.ReadNextBase();
        }

        return hostStrand;
    }

    /// <summary>A registration shift, performed blind to length: the probe is slid one base at a time until it pairs cleanly.</summary>
    public static bool FindsSeatingOn(ISequence probe, ISequence template)
    {
        bool seated;

        while (!(seated = probe.ClosesCleanlyAgainst(template)) && !probe.RunsPastEndOf(template))
            Shift(probe);

        return seated;
    }

    /// <summary>Transcribes from the CURRENT head position to the next 5hmu boundary (or the strand end)</summary>
    public static RNA TranscribeUntilBorder(ISequence template, IGenome genome)
    {
        RNA transcript = genome.CreateNewRnaSequence();

        while (template.HasNextBase() && template.Base != Base.Border)
        {
            transcript.Sequence.AddBase(GetWatsonAndCrickPair(template.Base));

            template.ReadNextBase();
        }
        return transcript;
    }

    /// <summary>A single shift: prepend one empty AP site, moving the strand one base downstream.</summary>
    public static void Shift(ISequence sequence)
    {
        sequence.Prepend(Base.Empty);
    }

    /// <summary>Releases a segment by stripping the leading boundary marks and shift overhang</summary>
    public static ISequence ReleaseSequence(ISequence sequence)
    {
        sequence.NewRead();
        while (sequence.Base == Base.Empty || sequence.Base == Base.Border)
            sequence.ReadNextBase();
        sequence.Cleave();

        return sequence;
    }

    /// <summary>The deterministic single-step isometry of the modulon vocabulary : A->C->G->T->A.</summary>
    public static Base AdvanceBase(Base current)
    {
        if (current == Base.A) return Base.C;
        if (current == Base.C) return Base.G;
        if (current == Base.G) return Base.T;
        if (current == Base.T) return Base.A;
        return current;
    }

    /// <summary>Materializes a stencil from its schema spec into an epigenetic guide (AP=accessible).</summary>
    public static ISequence Stencil(string stencilSpec, IGenome genome)
    {
        ISequence stencil = genome.CreateNewSequence();
        foreach (char mark in stencilSpec)
            stencil.AddBase(mark == 'e' ? Base.Empty : Base.Border);
        return stencil;
    }

    /// <summary>Advances the codeword where the target is accessible per stencil.</summary>
    public static ISequence AdvanceUnderStencil(ISequence target, ISequence stencil)
    {
        target.NewRead();
        stencil.NewRead();

        while (target.HasNextBase() && stencil.HasNextBase())
        {
            if (stencil.Base == Base.Empty) target.Base = AdvanceBase(target.Base);

            target.ReadNextBase();
            stencil.ReadNextBase();
        }

        return target;
    }

    /// <summary>Balances a pool to copies of each distinct type. Identity is Watson-Crick pairing.</summary>
    public static Pool<ISequence> Balance(Pool<ISequence> poolToBalance, int enrichment, IGenome genome)
    {
        Pool<Pool<ISequence>> typeGroups = new Pool<Pool<ISequence>>();

        while (poolToBalance.IsNotEmpty())
        {
            ISequence sequence = poolToBalance.PullOne();
            DoubleStrand sequenceDuplex = SynthesizeComplementaryStrand(sequence, genome);

            // find the type this sequence pairs with, draining then recycling so none are lost
            Pool<Pool<ISequence>> inspected = new Pool<Pool<ISequence>>();

            Pool<ISequence> matchingType = new Pool<ISequence>();

            while (typeGroups.IsNotEmpty())
            {
                Pool<ISequence> candidateType = typeGroups.PullOne();

                if (matchingType.IsEmpty())
                {
                    sequenceDuplex.NewRead();
                    if (sequenceDuplex.ClosesCleanlyAgainst(SynthesizeComplementaryStrand(candidateType.First(), genome)))
                        matchingType = candidateType;
                }
                inspected.Add(candidateType);
            }

            typeGroups.Recycle(inspected);

            if (matchingType.IsEmpty()) 
            { 
                matchingType = new Pool<ISequence>(enrichment); 
                typeGroups.Add(matchingType); 
            }

            matchingType.Add(sequence);
        }

        Pool<ISequence> balanced = new Pool<ISequence>();

        while (typeGroups.IsNotEmpty())
        {
            Pool<ISequence> typeGroup = typeGroups.PullOne();

            while (!typeGroup.Enriched()) typeGroup.Add(CopyStrand(typeGroup.First(), genome));

            balanced.Recycle(typeGroup);
        }

        return balanced.Shuffle(new Random(Settings.ShuffleSeed)); 
    }

    /// <summary>Enriches a pool by replicating each sequence times.</summary>
    public static Pool<ISequence> Enrich(Pool<ISequence> sequences, int rate, IGenome genome)
    {
        Pool<ISequence> enriched = new Pool<ISequence>();
        lock (sequences)
        {
            Pool<ISequence> originals = new Pool<ISequence>();
            while (sequences.IsNotEmpty())
            {
                ISequence sequence = sequences.PullOne();

                Pool<ISequence> copies = new Pool<ISequence>();
                while (copies.Count < rate) copies.Add(CopyStrand(sequence, genome));
                enriched.Recycle(copies);

                originals.Add(sequence);
            }
            sequences.Recycle(originals);   // restore the input pool
        }
        return enriched.Shuffle(new Random(Settings.ShuffleSeed));
    }
    
    /// <summary>Enriches a pool of duplexes by replicating both strands of each duplex times.</summary>
    public static Pool<DoubleStrand> Enrich(Pool<DoubleStrand> duplexes, int rate, IGenome genome)
    {
        Pool<DoubleStrand> enriched = new Pool<DoubleStrand>();
        lock (duplexes)
        {
            Pool<DoubleStrand> originals = new Pool<DoubleStrand>();
            while (duplexes.IsNotEmpty())
            {
                DoubleStrand duplex = duplexes.PullOne();

                Pool<DoubleStrand> copies = new Pool<DoubleStrand>();
                while (copies.Count < rate)
                    copies.Add(new DoubleStrand(CopyStrand(duplex.SenseStrand, genome),
                                                CopyStrand(duplex.ComplementaryStrand, genome)));
                enriched.Recycle(copies);

                originals.Add(duplex);
            }
            duplexes.Recycle(originals);   // restore the input pool
        }
        return enriched;
    }

}

