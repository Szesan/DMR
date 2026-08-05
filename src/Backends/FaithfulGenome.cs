using System;
using System.Collections.Generic;
using System.Linq;

// Genome medium as one FaithfulSequence, operated only through the reading head; spotted with free slots.
// Store occupies a free slot, fetch reads occupied frames, eliminate rewinds onto the mark and frees it.

public class FaithfulGenome : IGenome
{
    readonly FaithfulSequence medium;

    static Base OccupiedMark => Base.A;   // the concrete base substituted for the 'R' of Settings.SlotOccupied

    public FaithfulGenome(ISequence spottedMedium)
    {
        medium = new FaithfulSequence(spottedMedium.ToString());
    }

    public override string ToString() => medium.ToString();

    public ISequence CreateNewSequence() => new FaithfulSequence("");
    public DoubleStrand CreateNewDnaSequence() => new DoubleStrand(new FaithfulSequence(""), new FaithfulSequence(""));
    public RNA CreateNewRnaSequence() => new RNA(new FaithfulSequence(""));

    static bool IsEmptySite(Base site) => site == Base.Empty;
    static bool IsFlankMark(Base site) => site == Base.Border || site == Base.Empty;
    static bool IsRealBase(Base site)  => site != Base.Empty && site != Base.Border;   // a mark or payload base

    // ---- STORE (distributed): occupy the next free slot with each framed copy ----
    public void StoreDistributed(Pool<ISequence> sequences, string leftFlank, string rightFlank)
    {
        var left  = new FaithfulSequence(leftFlank);
        var right = new FaithfulSequence(rightFlank);
        var pending = new Pool<ISequence>(sequences);
        while (pending.IsNotEmpty())
            OccupyNextFreeSlot(pending.PullOne(), left, right);
    }

    public void StoreClustered(Pool<ISequence> sequences, string leftFlank, string rightFlank)
        => throw new NotSupportedException("FaithfulGenome does not store clustered.");

    void OccupyNextFreeSlot(ISequence payload, ISequence leftFlank, ISequence rightFlank)
    {
        medium.Direction = ReadingComplexDirection.Forward;
        medium.NewRead();
        while (medium.HasNextBase())
        {
            if (IsEmptySite(medium.Base) && HeadRestsOnFreeSlotMark())
            {
                medium.Base = OccupiedMark;                 // free -> occupied
                medium.ReadNextBase();
                LayDown(leftFlank);
                LayDown(payload);
                LayDown(rightFlank);
                return;
            }
            medium.ReadNextBase();
        }
    }

    // Walk the free motif alongside the empty run; head rests on the mark cell iff the whole motif matched.
    bool HeadRestsOnFreeSlotMark()
    {
        var motif = new FaithfulSequence(Settings.SlotFreeMotif);
        motif.NewRead();
        motif.ReadNextBase();                               // first empty already under the head
        while (motif.HasNextBase())
        {
            motif.ReadNextBase();
            medium.ReadNextBase();
            if (!medium.HasNextBase() || !IsEmptySite(medium.Base)) return false;
        }
        return true;
    }

    // lay a strand down at the head, advancing over the medium as each base is written
    void LayDown(ISequence strand)
    {
        strand.NewRead();
        while (strand.HasNextBase()) { medium.Base = strand.ReadNextBase(); medium.ReadNextBase(); }
    }

    public void NewRead()
    {
       
    }

    // ---- FETCH: read back each occupied frame's payload (mark + left flank ... right flank) ----
    public Pool<ISequence> FetchSequencesByFlankMotifs(string leftFlank, string rightFlank)
        => FetchFramed(new FaithfulSequence(leftFlank));

    public Pool<ISequence> FetchClustered(string leftFlank, string rightFlank)
        => FetchFramed(new FaithfulSequence(leftFlank));

    Pool<ISequence> FetchFramed(ISequence leftFlank)
    {
        var found = new Pool<ISequence>();
        medium.Direction = ReadingComplexDirection.Forward;
        medium.NewRead();
        while (medium.HasNextBase())
        {
            if (IsRealBase(medium.Base))                    // a candidate occupancy mark
            {
                medium.ReadNextBase();
                if (FlankConsumedAtHead(leftFlank))
                    found.Add(ReadPayloadToRightFlank());
            }
            else medium.ReadNextBase();
        }
        return found;
    }

    // Read to the right flank: stop at a mark-RUN; a lone internal border ('h' then content) stays in payload.
    FaithfulSequence ReadPayloadToRightFlank()
    {
        var payload = new FaithfulSequence("");
        while (medium.HasNextBase())
        {
            if (IsFlankMark(medium.Base))
            {
                var mark = medium.ReadNextBase();
                if (medium.HasNextBase() && IsFlankMark(medium.Base)) return payload;   // mark-run: the right flank
                payload.AddBase(mark);                                                  // lone internal border
            }
            else { payload.AddBase(medium.Base); medium.ReadNextBase(); }
        }
        return payload;
    }

    // consume a flank/motif from the head if it matches base-for-base, leaving the head just past it
    bool FlankConsumedAtHead(ISequence motif)
    {
        motif.NewRead();
        while (motif.HasNextBase())
        {
            if (!medium.HasNextBase() || medium.Base != motif.ReadNextBase()) return false;
            medium.ReadNextBase();
        }
        return true;
    }

    // ELIMINATE: match flank (+ optional payload prefix) forward, rewind onto the mark, flip it to free.
    public void EliminateByMotif(string leftFlankAndPrefix, string rightFlank)
    {
        var (leftFlank, payloadPrefix) = SplitFlankAndPrefix(leftFlankAndPrefix);

        medium.Direction = ReadingComplexDirection.Forward;
        medium.NewRead();
        while (medium.HasNextBase())
        {
            if (IsRealBase(medium.Base))
            {
                medium.ReadNextBase();
                if (FlankConsumedAtHead(leftFlank) && FlankConsumedAtHead(payloadPrefix))
                {
                    medium.Direction = ReadingComplexDirection.Reverse;
                    medium.ReadNextBase();                                              // step back off the payload
                    while (medium.HasNextBase() && IsRealBase(medium.Base)) medium.ReadNextBase();  // back over prefix
                    while (medium.HasNextBase() && IsFlankMark(medium.Base)) medium.ReadNextBase();  // back over flank
                    if (medium.HasNextBase() && IsRealBase(medium.Base)) medium.Base = Base.Empty;   // occupied -> free
                    medium.Direction = ReadingComplexDirection.Forward;
                    medium.ReadNextBase();
                }
            }
            else medium.ReadNextBase();
        }
    }

    // Split "flank + modulon" into (flank, prefix); a bare flank yields an empty prefix matching every frame.
    (ISequence flank, ISequence prefix) SplitFlankAndPrefix(string leftFlankAndPrefix)
    {
        foreach (var known in new[] { Settings.ModulonLeftFlank, Settings.RecipeLeftFlank, Settings.MapLeftFlank })
            if (leftFlankAndPrefix.StartsWith(known, StringComparison.Ordinal))
                return (new FaithfulSequence(known),
                        new FaithfulSequence(leftFlankAndPrefix.Substring(known.Length)));
        return (new FaithfulSequence(leftFlankAndPrefix), new FaithfulSequence(""));
    }

    // Damage: environment mutates a site; walk the head there (no indexing), sparing marks.
    public void ApplyMutation(int position, Base mutatedBase)
    {
        medium.Direction = ReadingComplexDirection.Forward;
        medium.NewRead();
        var steps = position;
        while (steps > 0 && medium.HasNextBase()) { medium.ReadNextBase(); steps--; }
        if (!medium.HasNextBase()) return;
        if (medium.Base == Base.Border) return;
        medium.Base = mutatedBase;
    }

    public Pool<TOut> ProcessAtScale<TIn, TOut>(Func<TIn, TOut> job, Pool<TIn> units)
    {
        var results = new Pool<TOut>();
        var pending = new Pool<TIn>(units);
        while (pending.IsNotEmpty()) results.Add(job(pending.PullOne()));
        return results;
    }
}
