// Damage scheduling, scoring and reporting. Harness, not biology. See README §3.
// Orchestrator: drives one preservation model under a seeded clustered-damage stream (see DAMAGE_ORCHESTRATOR_PLAN.md).

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Diagnostics; 

public enum PreservationModelKind { FullDmr, Dmr, Naive }

public class Simulation
{
    readonly PreservationModelKind kind;
    readonly int seed;
    readonly string outputDirectory;

    GenomeFactory factory;
    IGenome genome;
    IPreservationModel model;
    HashSet<string> ancestralGenes;
    Random damageRng;
    int ticksApplied;

    public Simulation(PreservationModelKind kind = PreservationModelKind.FullDmr, int seed = 100,
                      string outputDirectory = "runs")
    {
        this.kind = kind;
        this.seed = seed;
        this.outputDirectory = outputDirectory;
    }

    string RunName => $"{kind}_seed{seed}_maint{Settings.MaintenanceInterval}_report{Settings.ReportInterval}";
    string ResultsPath => Path.Combine(outputDirectory, $"{RunName}.csv");
    string GenomePath  => Path.Combine(outputDirectory, $"{RunName}.genome.txt");
    string AncestralPath => Path.Combine(outputDirectory, $"{RunName}.ancestral.txt");

    // ---- construction ----

    public void Build()
    {
        factory = new GenomeFactory(seed);
        ancestralGenes = new HashSet<string>(factory.AncestralGenes.Select(gene => gene.ToString()));
        genome = NewGenome();
        model = NewModel(genome);
        damageRng = new Random(seed);
        ticksApplied = 0;

        DumpAncestralReference();
    }

    IGenome NewGenome() => kind == PreservationModelKind.Naive
        ? factory.BuildNaiveGenome(new PackedGenome(factory.NewBlankMedium()))
        : factory.BuildDmrGenome(new PackedGenome(factory.NewBlankMedium()));

    IPreservationModel NewModel(IGenome builtGenome) => kind switch
    {
        PreservationModelKind.FullDmr => new Dmr(builtGenome),
        PreservationModelKind.Dmr     => new DmrNoMaintenance(builtGenome),
        _                             => new NaiveDuplication(builtGenome),
    };

    // ---- clustered damage: one hit rolled once, applied to this run's genome ----

    static readonly Base[] Alphabet = { Base.A, Base.C, Base.G, Base.T };

    void Hit(int length, int position)
    {
        if (position < 0 || position >= length) return;      // out of bounds consumes no roll
        var replacement = Alphabet[damageRng.Next(Alphabet.Length)];
        genome.ApplyMutation(position, replacement);         // each backend spares its own epigenetic marks
    }

    void DamageEvent(int length)
    {
        if (damageRng.NextDouble() >= Settings.ClusterInitiationProbability) return;
        var centre = damageRng.Next(length);
        Hit(length, centre);
        var offset = -Settings.ScatterRadius;
        while (offset <= Settings.ScatterRadius)
        {
            if (offset != 0 && damageRng.NextDouble() < Settings.ScatterProbability) Hit(length, centre + offset);
            offset++;
        }
    }

    // Advance the RNG without touching a genome, so a resumed run continues the identical stream.
    void ReplayDamage(int ticks, int length)
    {
        var replayed = 0;
        while (replayed < ticks)
        {
            if (damageRng.NextDouble() < Settings.ClusterInitiationProbability)
            {
                var centre = damageRng.Next(length);
                damageRng.Next(Alphabet.Length);                 // centre is always in bounds
                var offset = -Settings.ScatterRadius;
                while (offset <= Settings.ScatterRadius)
                {
                    if (offset != 0 && damageRng.NextDouble() < Settings.ScatterProbability)
                    {
                        var position = centre + offset;
                        if (position >= 0 && position < length) damageRng.Next(Alphabet.Length);
                    }
                    offset++;
                }
            }
            replayed++;
        }
    }

    // ---- measurement ----

    (int recovered, int spurious) Score(List<ISequence> presented)
    {
        var present = new HashSet<string>(presented.Select(gene => gene.ToString()));
        var recovered = ancestralGenes.Count(gene => present.Contains(gene));
        var spurious = present.Count(gene => !ancestralGenes.Contains(gene));
        return (recovered, spurious);
    }

    // ---- the run ----

    public void Run(int ticks, double recoveryFloor = 0.0)
    {
        if (model == null) Build();
        var length = genome.ToString().Length;
        WriteHeaderIfMissing();

        var target = ticksApplied + ticks;
        var collapsed = false;

        var wholeRun = Stopwatch.StartNew();
        var reportedCycle = new Stopwatch();

        while (ticksApplied < target && !collapsed)
        {
            ticksApplied++;
            DamageEvent(length);

            reportedCycle.Restart();
            var presented = model.RunCycle(ticksApplied);
            reportedCycle.Stop();

            // Checkpoint after each refresh so a long run interrupted mid-flight stays resumable.
            if (Settings.CheckpointOnMaintenance && ticksApplied.DueForMaintenance()) DumpGenome();

            if (ticksApplied.DueForReport())
            {
                var (recovered, spurious) = Score(presented);
                var cycleMs = reportedCycle.ElapsedMilliseconds;

                AppendRow(ticksApplied, recovered, spurious, cycleMs);
                Console.WriteLine($"  tick {ticksApplied,6}  recovered {recovered,3}/{ancestralGenes.Count}" +
                                  $"  spurious {spurious}  cycle {cycleMs} ms");

                collapsed = (double)recovered / ancestralGenes.Count < recoveryFloor;
            }
        }

        wholeRun.Stop();

        var reason = collapsed ? $"recovery fell below {recoveryFloor:P0}" : "reached tick budget";
        File.AppendAllText(ResultsPath,
            $"# ended at tick {ticksApplied}: {reason}; run took {wholeRun.ElapsedMilliseconds} ms\n");
        Console.WriteLine($"  ended at tick {ticksApplied}: {reason}; run took {wholeRun.ElapsedMilliseconds} ms");

        DumpGenome();
    }

    void WriteHeaderIfMissing()
    {
        if (File.Exists(ResultsPath)) return;
        Directory.CreateDirectory(outputDirectory);
        File.WriteAllText(ResultsPath,
            $"# model={kind} seed={seed} maintenanceInterval={Settings.MaintenanceInterval} reportInterval={Settings.ReportInterval}\n" +
            $"# genomeSize={Settings.GenomeSize} genes={Settings.GeneCount} exonLength={Settings.ExonLength} exonsPerGene={Settings.ExonsPerGene}\n" +
            "tick,recovered,spurious,totalGenes,cycleMs\n");
    }

    void AppendRow(int tick, int recovered, int spurious, long cycleMs)
        => File.AppendAllText(ResultsPath, $"{tick},{recovered},{spurious},{ancestralGenes.Count},{cycleMs}\n");

    // Pristine genes and exon bodies, written before any damage: the independent reference.
    void DumpAncestralReference()
    {
        Directory.CreateDirectory(outputDirectory);
        var reference = new StringBuilder();
        reference.AppendLine($"# model={kind} seed={seed} genes={Settings.GeneCount} " +
                             $"exonsPerGene={Settings.ExonsPerGene} exonLength={Settings.ExonLength} " +
                             $"modulonLength={Settings.ModulonLength}");

        reference.AppendLine("# exon types: modulon<TAB>body");
        var modulons = factory.CanonicalModulons;
        var exonTypes = factory.CanonicalExonBodies;
        for (int type = 0; type < exonTypes.Count; type++)
            reference.AppendLine($"exon\t{modulons[type]}\t{exonTypes[type]}");

        reference.AppendLine("# genes: index<TAB>gene<TAB>exon bodies in order");
        var geneExons = factory.AncestralGeneExons;
        for (int gene = 0; gene < factory.AncestralGenes.Count; gene++)
            reference.AppendLine($"gene\t{gene}\t{factory.AncestralGenes[gene]}\t" +
                                 string.Join(" ", geneExons[gene].Select(exon => exon.ToString())));

        File.WriteAllText(AncestralPath, reference.ToString());
    }

    // Genome plus tick count is the whole resumable state; temp-then-move so a partial write can't truncate it.
    void DumpGenome()
    {
        var pending = GenomePath + ".tmp";
        File.WriteAllText(pending, $"{ticksApplied}\n{genome}\n");
        File.Move(pending, GenomePath, overwrite: true);
    }

    // Reload a dumped genome and fast-forward the RNG so the stream continues unbroken.
    public void Resume()
    {
        Build();
        if (!File.Exists(GenomePath)) return;

        var dump = File.ReadAllLines(GenomePath);
        ticksApplied = int.Parse(dump[0]);
        genome = new PackedGenome(new PackedSequence(dump[1]));
        model = NewModel(genome);
        damageRng = new Random(seed);
        ReplayDamage(ticksApplied, genome.ToString().Length);
    }
}