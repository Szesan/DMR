public class DoubleStrand : ISequence
{   
    /// <summary>Gets the primary or sense strand of the DNA molecule.</summary>
    public ISequence SenseStrand { get; }

    /// <summary>Gets the complementary or anti-sense strand of the DNA molecule.</summary>
    public ISequence ComplementaryStrand { get; }


    public DoubleStrand(ISequence senseStrand, ISequence complementaryStrand)
    {
        SenseStrand = senseStrand;
        ComplementaryStrand = complementaryStrand;
    }

    // Direction propagates to both strands so the duplex reads in lockstep; forward is relative to sense.
    public ReadingComplexDirection Direction
    {
        get => SenseStrand.Direction;
        set { SenseStrand.Direction = value; ComplementaryStrand.Direction = value; }
    }

    public void NewRead()
    {
        SenseStrand.NewRead();
        ComplementaryStrand.NewRead();
    }

    public bool HasNextBase()
    {
        return SenseStrand.HasNextBase() && ComplementaryStrand.HasNextBase();
    }

    public bool RunsPastEndOf(ISequence template)
    {
        return SenseStrand.RunsPastEndOf(template);
    }

    public bool ContainsNoneOf(Base target)
    {
        return SenseStrand.ContainsNoneOf(target);
    }

    public Base Base
    {
        get => SenseStrand.Base;
        set
        {
            SenseStrand.Base = value;
            ComplementaryStrand.Base = Enzyme.GetWatsonAndCrickPair(value);
        }
    }

    public Base ReadNextBase()
    {
        ComplementaryStrand.ReadNextBase();
        return SenseStrand.ReadNextBase();
    }

    public void AddBase(Base baseToAdd)
    {
        SenseStrand.AddBase(baseToAdd);
        ComplementaryStrand.AddBase(Enzyme.GetWatsonAndCrickPair(baseToAdd));        
    }

    public void Prepend(Base baseToPrepend)
    {
        SenseStrand.Prepend(baseToPrepend);
        ComplementaryStrand.Prepend(Enzyme.GetWatsonAndCrickPair(baseToPrepend)); 
    }

    public ISequence CreateEmptyBackbone()
    {
        return new DoubleStrand(
            SenseStrand.CreateEmptyBackbone(),
            ComplementaryStrand.CreateEmptyBackbone()
        );
    }

    public ISequence Cleave()
    {
        return new DoubleStrand(
            SenseStrand.Cleave(),
            ComplementaryStrand.Cleave()
        );
    }

    public ISequence Ligate(ISequence other)
    {
        if (other is DoubleStrand dsOther)
        {
            SenseStrand.Ligate(dsOther.SenseStrand);
            ComplementaryStrand.Ligate(dsOther.ComplementaryStrand);
        }
        else
        {
            SenseStrand.Ligate(other);           
            ComplementaryStrand.Ligate(other.CreateEmptyBackbone());
        }

        return this;
    }

    // Closure for a DoubleStrand. Biologically, two duplexes anneal if BOTH strand pairings hold:
    public bool ClosesCleanlyAgainst(ISequence template)
    {
        if (template is DoubleStrand other)
            return SenseStrand.ClosesCleanlyAgainst(other.ComplementaryStrand)
                && ComplementaryStrand.ClosesCleanlyAgainst(other.SenseStrand);

        return SenseStrand.ClosesCleanlyAgainst(template);
    }

    // A decoder column closes here if: there is a base on both heads to close against; 
    public bool ClosesAtHead(DoubleStrand other)
    {
        if (!HasNextBase() || !other.HasNextBase()) return false;
        if (!Enzyme.BasePairsWith(SenseStrand.Base, other.ComplementaryStrand.Base)) return false;
        if (SenseStrand.Base != Base.Empty && other.SenseStrand.Base != Base.Empty) return false;
        return true;
    }

    public override string ToString()
    {
        return SenseStrand.ToString();
    }
}

public class RNA : ISequence
{
    public ISequence Sequence { get; }

    public RNA(ISequence sequence)
    {
        Sequence = sequence;
    }

    // Every operation delegates to the transcribed strand; RNA is a single-stranded carrier.
    public ReadingComplexDirection Direction
    {
        get => Sequence.Direction;
        set => Sequence.Direction = value;
    }

    public Base Base
    {
        get => Sequence.Base;
        set => Sequence.Base = value;
    }

    public void NewRead()
    {
        Sequence.NewRead();
    }

    public bool HasNextBase()
    {
        return Sequence.HasNextBase();
    }

    public bool RunsPastEndOf(ISequence template)
    {
        return Sequence.RunsPastEndOf(template);
    }

    public bool ContainsNoneOf(Base target)
    {
        return Sequence.ContainsNoneOf(target);
    }

    public Base ReadNextBase()
    {
        return Sequence.ReadNextBase();
    }

    public void AddBase(Base baseToAdd)
    {
        Sequence.AddBase(baseToAdd);
    }

    public void Prepend(Base baseToPrepend)
    {
        Sequence.Prepend(baseToPrepend);
    }

    public ISequence CreateEmptyBackbone()
    {
        return Sequence.CreateEmptyBackbone();
    }

    public ISequence Cleave()
    {
        return Sequence.Cleave();
    }

    public ISequence Ligate(ISequence other)
    {
        return Sequence.Ligate(other);
    }

    public bool ClosesCleanlyAgainst(ISequence template)
    {
        return Sequence.ClosesCleanlyAgainst(template);
    }

    public override string ToString()
    {
        return Sequence.ToString();
    }
}
