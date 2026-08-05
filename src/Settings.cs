// All configuration. See README §11.
public static class Settings
{
    public const char Terminator = '\u0000';   // sensed past a strand's end 
    public const int  Start = 0;                // the reading head's origin

    // ---- Motif widths ----
    public static readonly int ModulonLength = 9;
    public static readonly int ArchiveWordLength = 8;
    // Exon body size; gene length = ExonsPerGene x ExonLength. Settable so benchmarks can sweep it.
    public static int ExonLength = 90;
    public static readonly int FirstPosition = 0;

    // ---- Enrichment (copies per stored type) ----
    public static readonly int RootEnrichmentRate = 12;        // candidates drawn at each gene position
    public static readonly int IdentifierEnrichmentRate = 12;  // copies per map / recipe
    public static readonly string ModulonLeftFlank   = "hhhhhhh"; // 7  roots' left flank
    public static readonly string RecipeLeftFlank     = "hhhhhhe"; // 6
    public static readonly string MapLeftFlank        = "hhhhhee"; // 5
    public static readonly string MapRightFlank       = "hhhheee"; // 4
    public static readonly string RecipeRightFlank    = "hhheeee"; // 3
    public static readonly string ModulonRightFlank   = "hheeeee"; // 2  roots' right flank

    public static readonly string SlotFreeMotif = "eeeeeee";

    public static readonly string SlotOccupied = "eeeeeeR"; 

    public static readonly string ModulonStencil = "hhheeehhh";               // 9:  shoulder | centre | shoulder
    public static readonly string MapStencil     = "hhhhhhhhhhheeehhh";        // 17: 11 protected (word + modulon prefix) | centre | shoulder

    // Genome-side scan span for one frame; tracks the ExonLength knob.
    public static int MaxPayloadLength => ModulonLength + ExonLength + 8;

    public static readonly int ExonTypeCount = 6;
    public static readonly int GeneCount = 16;
    public static readonly int ExonsPerGene = 5;

    public static readonly int EnrichmentAbundant = 30;

    public static readonly int NaiveCopiesPerGene = 3;

    public static readonly int ScatterStride = 200;

    public static int GenomeSize = 60000;

    // Maintenance cadence: swap-and-amplify runs when time.DueForMaintenance().
    public static readonly int MaintenanceInterval = 600;

    public static int ReportInterval = 50;
    // Checkpoint the medium at each refresh so an interrupted run stays resumable.
    public static bool CheckpointOnMaintenance = true;
    public static double ClusterInitiationProbability = 0.25;
    public static double ScatterProbability = 0.50;
    public static int ScatterRadius = 15;
    public static readonly int ShuffleSeed = 12345;

    
    public static readonly int ModulonMinSeparation = 4;
    public static readonly int ArchiveWordMinSeparation = 3;

    // ---- Diagnostics ----
    public static readonly string DebugLogPath = "logs/DebugLog.txt";
    public static bool DebugLogEnabled = true;
}
