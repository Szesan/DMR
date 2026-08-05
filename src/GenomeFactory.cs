// Builds the initial state: mints one ancestral gene set from reusable exon parts, then
// stores it as a scattered DMR archive or as clustered whole genes. Construction edge, see README §3.
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text; 

    // Factory knobs. Structural constants (flank motifs) are not knobs and stay on Settings.
public class FactoryKnobs
{
    public int RootEnrichmentRate;
    public int IdentifierEnrichmentRate;
    public int ExonTypeCount;
    public int ExonsPerGene;
    public int GeneCount;
    public int ExonLength;
    public int ModulonLength;
    public int ArchiveWordLength;
    public int NaiveCopiesPerGene;
    public int GenomeSize;

    public static FactoryKnobs FromParams() => new FactoryKnobs
    {
        RootEnrichmentRate       = Settings.RootEnrichmentRate,
        IdentifierEnrichmentRate = Settings.IdentifierEnrichmentRate,
        ExonTypeCount            = Settings.ExonTypeCount,
        ExonsPerGene             = Settings.ExonsPerGene,
        GeneCount                = Settings.GeneCount,
        ExonLength               = Settings.ExonLength,
        ModulonLength            = Settings.ModulonLength,
        ArchiveWordLength        = Settings.ArchiveWordLength,
        NaiveCopiesPerGene       = Settings.NaiveCopiesPerGene,
        GenomeSize               = Settings.GenomeSize,
    };
}

public class GenomeFactory
{
    readonly char[] Bases = { 'A', 'C', 'G', 'T' };
    const char Internal = 'h';   // the root's internal modulon|exon boundary (payload, not flank)

    static readonly string[] Modulons =
    {
        "CTTGTTCCC", "CTCTTGTTC", "GAGATTTCA", "GAGAGGAAC", "AGAGATTTC", "AAATCTCTA"
    };
    static readonly string[] ArchiveWords =
    {
        "ACACTAAG", "AAAAGCTT", "AATTGACC", "AACATGAT", "ATTTGGAG", "AAACGTGA"
    };

    readonly Random rng;

    // Tunable knobs captured at construction — from the Settings file by default, or overridden for tests.
    readonly FactoryKnobs P;

    // minted vocabularies / structure, indexed by exon TYPE (0..ExonTypeCount-1)
    readonly List<string> exons = new List<string>();        // exon bodies (gene-roots)
    readonly List<string> modulons = new List<string>();     // operational tags
    readonly List<string> archiveWords = new List<string>(); // stable storage ids
    readonly List<List<int>> geneTypeOrders = new List<List<int>>(); // per gene: type sequence

    // the held ancestral reference: assembled gene sequences (exon bodies concatenated)
    readonly List<ISequence> ancestralGenes = new List<ISequence>();

    // Default: read all knobs from the Settings file.
    public GenomeFactory(int seed) : this(seed, FactoryKnobs.FromParams()) { }

    // Test override: take the knobs from the passed config instead of the Settings file.
    public GenomeFactory(int seed, FactoryKnobs knobs)
    {
        P = knobs;
        rng = new Random(seed);
        MintExons();
        MintIdentifierVocabularies();
        BuildGenes();
    }

    // The shared reference the orchestrator scores every model against.
    public IReadOnlyList<ISequence> AncestralGenes => ancestralGenes;

    // Canonical distinct payloads per layer, pre-enrichment. Exposed for scoring.
    public IReadOnlyList<string> CanonicalRootPayloads   => RootPayloads();    // 1 per exon type
    public IReadOnlyList<string> CanonicalMapPayloads    => MapPayloads();     // 1 per exon type
    public IReadOnlyList<string> CanonicalRecipePayloads => RecipePayloads();  // 1 per gene
    public IReadOnlyList<string> CanonicalModulons       => modulons;          // the seating key per exon type
    public IReadOnlyList<string> CanonicalExonBodies     => exons;             // 1 body per exon type
    public int ModulonWidth => P.ModulonLength;

    // Each ancestral gene's exon bodies in assembly order; bodies are the primary minted objects.
    public IReadOnlyList<IReadOnlyList<ISequence>> AncestralGeneExons
    {
        get
        {
            var perGene = new List<IReadOnlyList<ISequence>>();
            foreach (var typeOrder in geneTypeOrders)
            {
                var bodies = new List<ISequence>();
                foreach (var type in typeOrder) bodies.Add(new FaithfulSequence(exons[type]));
                perGene.Add(bodies);
            }
            return perGene;
        }
    }

    // The shared medium size — identical across models (equal budget).
    public int GenomeSize => P.GenomeSize;

    // Build the scattered DMR genome from the ancestral decomposition.

    // Backend-swappable: populates any blank IGenome with the same layout.
    public IGenome BuildDmrGenome(IGenome genome)
    {
        StoreRoots(genome);
        StoreMaps(genome);
        StoreRecipes(genome);

        return genome;
    }

    // Build the naive genome: whole genes stored clustered, so copies share a fate under correlated damage.
    public IGenome BuildNaiveGenome(IGenome genome)
    {
        var copies = new List<ISequence>();
        foreach (var gene in ancestralGenes)
            for (var k = 0; k < P.NaiveCopiesPerGene; k++)
                copies.Add(new FaithfulSequence(gene.ToString()));

        genome.StoreClustered(new Pool<ISequence>(copies), Settings.RecipeLeftFlank, Settings.RecipeRightFlank);
        return genome;
    }

    static void Shuffle(List<ISequence> xs, Random rng)
    {
        for (int i = xs.Count - 1; i > 0; i--) { int j = rng.Next(i + 1); (xs[i], xs[j]) = (xs[j], xs[i]); }
    }

    // The distinct payloads of each layer, in type/gene order (schema made legible).
    List<string> RootPayloads() =>
        Enumerable.Range(0, P.ExonTypeCount).Select(t => modulons[t] + Internal + exons[t]).ToList();
    List<string> MapPayloads() =>
        Enumerable.Range(0, P.ExonTypeCount).Select(t => archiveWords[t] + modulons[t]).ToList();
    List<string> RecipePayloads() =>
        geneTypeOrders.Select(order => string.Concat(order.Select(t => archiveWords[t]))).ToList();

    // A fresh blank medium (GenomeSize random bases) the caller can wrap in any backend genome.
    public ISequence NewBlankMedium() => BlankMedium();

    // Minting.
    void MintExons()
    {
        var seen = new HashSet<string>();
        while (exons.Count < P.ExonTypeCount)
        {
            var body = RandomBases(P.ExonLength);
            if (seen.Add(body)) exons.Add(body);
        }
    }

    // Identifiers come from the curated comma-free table in Settings, not random draws: random
    // codewords collide across the frame seam.
    void MintIdentifierVocabularies()
    {
        DrawVocabulary(modulons, Modulons, P.ModulonLength, "modulon");
        DrawVocabulary(archiveWords, ArchiveWords, P.ArchiveWordLength, "archive word");
    }

    // Copy a codeword table into the vocabulary, asserting it fits the schema. Fails loudly here.
    void DrawVocabulary(List<string> into, string[] table, int length, string what)
    {
        if (table.Length < P.ExonTypeCount)
            throw new InvalidOperationException(
                $"Curated {what} table has {table.Length} codewords; need {P.ExonTypeCount}.");
        if (table.Take(P.ExonTypeCount).Distinct().Count() != P.ExonTypeCount)
            throw new InvalidOperationException($"Curated {what} table has duplicate codewords.");
        if (table.Take(P.ExonTypeCount).Any(w => w.Length != length))
            throw new InvalidOperationException($"Curated {what} table has a codeword of wrong length (expected {length}).");

        for (var t = 0; t < P.ExonTypeCount; t++) into.Add(table[t]);
    }

    void BuildGenes()
    {
        if (P.ExonsPerGene > P.ExonTypeCount)
        {
            throw new InvalidOperationException("Cannot build genes with unique exons: ExonsPerGene exceeds ExonTypeCount.");
        }

    // A gene is its ordered exon types, so GeneCount must not exceed P(ExonTypeCount, ExonsPerGene).
        if (P.GeneCount > DistinctOrderingCount(P.ExonTypeCount, P.ExonsPerGene))
            throw new InvalidOperationException(
                $"Cannot build {P.GeneCount} distinct genes: only {DistinctOrderingCount(P.ExonTypeCount, P.ExonsPerGene)} " +
                $"orderings exist for {P.ExonsPerGene} of {P.ExonTypeCount} types.");

        for (var g = 0; g < P.GeneCount; g++)
        {
            List<int> typeOrder;

            do
            {
                typeOrder = new List<int>();
                for (var slot = 0; slot < P.ExonsPerGene; slot++)
                {
                    int type;
                    do { type = rng.Next(P.ExonTypeCount); } while (typeOrder.Contains(type));
                    typeOrder.Add(type);
                }
            } while (geneTypeOrders.Any(existing => existing.SequenceEqual(typeOrder)));  // no duplicate gene

            var assembled = new StringBuilder();
            foreach (var type in typeOrder) assembled.Append(exons[type]);

            geneTypeOrders.Add(typeOrder);
            ancestralGenes.Add(new FaithfulSequence(assembled.ToString()));
        }
    }

    // Number of distinct ordered exon-type sequences: P(types, perGene) = types! / (types - perGene)!
    static int DistinctOrderingCount(int types, int perGene)
    {
        var count = 1;
        for (var k = 0; k < perGene; k++) count *= (types - k);
        return count;
    }

    // Build payload -> balance -> store. The genome frames on store.
    void StoreRoots(IGenome genome)
    {
        // one root payload per TYPE: [modulon] h [exon]; scattered under roots' flanks (7 .. 2).
        var roots = new List<ISequence>();
        for (var t = 0; t < P.ExonTypeCount; t++)
            roots.Add(new FaithfulSequence(modulons[t] + Internal + exons[t]));
        genome.StoreDistributed(Enzyme.Balance(new Pool<ISequence>(roots), P.RootEnrichmentRate, genome),
                                Settings.ModulonLeftFlank, Settings.ModulonRightFlank);
    }

    void StoreMaps(IGenome genome)
    {
        // one map payload per TYPE: [archiveWord][modulon]; scattered under map flanks (5 .. 4).
        var maps = new List<ISequence>();
        for (var t = 0; t < P.ExonTypeCount; t++)
            maps.Add(new FaithfulSequence(archiveWords[t] + modulons[t]));
        genome.StoreDistributed(Enzyme.Balance(new Pool<ISequence>(maps), P.IdentifierEnrichmentRate, genome),
                                Settings.MapLeftFlank, Settings.MapRightFlank);
    }

    void StoreRecipes(IGenome genome)
    {
        // one recipe payload per GENE: [archiveWord]...; scattered under recipe flanks (6 .. 3).
        var recipes = new List<ISequence>();
        foreach (var typeOrder in geneTypeOrders)
        {
            var payload = new StringBuilder();
            foreach (var type in typeOrder)
                payload.Append(archiveWords[type]);
            recipes.Add(new FaithfulSequence(payload.ToString()));
        }
        genome.StoreDistributed(Enzyme.Balance(new Pool<ISequence>(recipes), P.IdentifierEnrichmentRate, genome),
                                Settings.RecipeLeftFlank, Settings.RecipeRightFlank);
    }

    // Helpers.
    ISequence BlankMedium()
    {
        var sb = new StringBuilder(P.GenomeSize);
        for (var i = 0; i < P.GenomeSize; i++) sb.Append(Bases[rng.Next(Bases.Length)]);
        return new FaithfulSequence(sb.ToString());
    }

    string RandomBases(int n)
    {
        var sb = new StringBuilder(n);
        for (var i = 0; i < n; i++) sb.Append(Bases[rng.Next(Bases.Length)]);
        return sb.ToString();
    }
}