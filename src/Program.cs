// Entry point: dotnet run -- [naive|dmr|fulldmr] [seed] [ticks] [floor] [--resume]

using System;
using System.Linq;

public static class Program
{
    static PreservationModelKind ParseModel(string name) => name?.ToLowerInvariant() switch
    {
        "naive"   => PreservationModelKind.Naive,
        "dmr"     => PreservationModelKind.Dmr,
        "fulldmr" => PreservationModelKind.FullDmr,
        null      => PreservationModelKind.FullDmr,
        _         => throw new ArgumentException($"unknown model '{name}' (expected naive, dmr or fulldmr)"),
    };

    public static void Main(string[] args)
    {
        var arguments = args.Where(argument => !argument.StartsWith("--")).ToList();
        var resuming  = args.Any(argument => argument == "--resume");

        var chosenModel = ParseModel(arguments.ElementAtOrDefault(0));
        var chosenSeed  = int.Parse(arguments.ElementAtOrDefault(1) ?? "100");
        var chosenTicks = int.Parse(arguments.ElementAtOrDefault(2) ?? "1000");
        var chosenFloor = double.Parse(arguments.ElementAtOrDefault(3) ?? "0.5");

        Console.WriteLine($"model={chosenModel} seed={chosenSeed} ticks={chosenTicks} " +
                          $"maintenance={Settings.MaintenanceInterval} report={Settings.ReportInterval} " +
                          $"floor={chosenFloor:P0} resume={resuming}");

        var simulation = new Simulation(chosenModel, chosenSeed);
        if (resuming) simulation.Resume(); else simulation.Build();
        simulation.Run(chosenTicks, chosenFloor);
    }
}
