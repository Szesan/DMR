// Pool<T>, logging and helpers. Harness plumbing.
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using static Settings;

public static class Logger
{
    private static readonly BlockingCollection<string> _logQueue = new BlockingCollection<string>();

    static Logger()
    {
        Task.Factory.StartNew(() =>
        {
            foreach (var message in _logQueue.GetConsumingEnumerable())
            {
                try
                {
                    var dir = Path.GetDirectoryName(Settings.DebugLogPath);
                    if (!string.IsNullOrEmpty(dir)) Directory.CreateDirectory(dir);
                    File.AppendAllText(Settings.DebugLogPath, message + Environment.NewLine);
                }
                catch { }
            }
        }, TaskCreationOptions.LongRunning);
    }

    public static void DebugLog(string message) 
    {
        if (!Settings.DebugLogEnabled) return;
        var stamp = DateTime.Now.ToString("yyyy-MM-dd HH:mm:ss.fff");
        _logQueue.Add($"[{stamp}] {message}");
    }
}

public class Pool<T> : Collection<T>
{
    // caps the population at EnrichmentRate (0 = unlimited)
    
    public int EnrichmentRate { get; set; }
 
    public Pool(int enrichmentRate = 0) : base() { EnrichmentRate = enrichmentRate; }
 
    public Pool(IEnumerable<T> items, int enrichmentRate = 0) : this(enrichmentRate)
    {
        foreach (var item in items) Add(item);
    }
 
    protected override void InsertItem(int index, T item)
    {
        if (EnrichmentRate == 0 || Count < EnrichmentRate) base.InsertItem(index, item);
    }
 
    public T PullOne() { var pulled = this.First(); Remove(pulled); return pulled; }
 
    public bool IsNotEmpty() => Count > 0;
 
    public bool IsEmpty() => Count == 0;
 
    public Pool<T> Shuffle(Random rng)
    {
        if (rng == null) throw new ArgumentNullException(nameof(rng));
        int n = Count;
        while (n > 1) { n--; int k = rng.Next(n + 1); (this[k], this[n]) = (this[n], this[k]); }
        return this;
    }
 
    // Re-orderings of this pool sharing members, not copies: callers needing to mutate copy first.
    public Pool<Pool<T>> ShuffledOrderings(int count, Random rng)
    {
        var orderings = new Pool<Pool<T>>();
        while (orderings.Count < count) orderings.Add(new Pool<T>(this).Shuffle(rng));
        return orderings;
    }
 
    public void Recycle(Pool<T> poolToRecycle)
    {
        if (poolToRecycle == null) return;
        while (poolToRecycle.IsNotEmpty()) Add(poolToRecycle.PullOne());
    }
 
    public List<T> ToList() => new List<T>(this);
 
    public bool Enriched() => Count == EnrichmentRate;
}

// Extension Methods

public static class Extensions
{

    public static bool DueForMaintenance(this int time)
    {
        if (time == 0) return false;

        return time % Settings.MaintenanceInterval == 0;
    }

    public static bool DueForReport(this int time)
    {
        if (time == 0) return false;

        return time % Settings.ReportInterval == 0;
    }

    public static void Reset(this ref int number)
    {
        number = 0;
    }

    public static void Move(this ref int number)
    {
        number++;
    }

    public static void MoveReverse(this ref int number)
    {
        number--;
    }

    public static char SymbolAtHead(this string strand, int head)
        => (uint)head < (uint)strand.Length ? strand[head] : Terminator;

    public static string AppendSymbol(this string strand, char symbol)  => strand + symbol;

    public static string PrependSymbol(this string strand, char symbol) => symbol + strand;

    public static string WriteSymbolAtHead(this string strand, int head, char symbol)
        => strand.Substring(0, head) + symbol + strand.Substring(head + 1);

    public static (string up, string down) SplitAtHead(this string strand, int head)
        => (strand.Substring(0, head), strand.Substring(head));
}
