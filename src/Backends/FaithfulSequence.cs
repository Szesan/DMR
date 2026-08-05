// Strict reading-head medium: the demonstration that the operations are physically possible.
// See README §4.1.
public partial class FaithfulSequence : ISimSequence
{
    const char  Adenine = 'A', 
                Cytosine = 'C', 
                Thymine = 'T', 
                Guanine = 'G', 
                Border = 'h', 
                Empty = 'e';

    string strand;

    int head;   // where the reading head sits; steps forward or backward per Direction, or resets

    public ReadingComplexDirection Direction { get; set; } = ReadingComplexDirection.Forward;

    public FaithfulSequence(string bases) 
    { 
        strand = bases ?? ""; 
        head = Settings.Start; 
    }

    char ToChar(Base b) => b switch
    {
        Base.A => Adenine, 
        Base.C => Cytosine, 
        Base.G => Guanine, 
        Base.T => Thymine,
        Base.Border => Border, 
        _ => Empty,
    };

    Base ToBase(char c) => c switch
    {
        Adenine => Base.A, Cytosine => Base.C, Guanine => Base.G, Thymine => Base.T,
        Border => Base.Border, _ => Base.Empty,
    };

    // ---- strict reading head ----
    public void NewRead() => head.Reset();   // origin-only; Direction does not change where a fresh read starts

    // Off the strand is sensed in the direction of travel.
    public bool HasNextBase() =>
        Direction == ReadingComplexDirection.Forward
            ? strand.SymbolAtHead(head) != Settings.Terminator
            : head >= Settings.Start && strand.SymbolAtHead(head) != Settings.Terminator;

    public Base Base
    {
        get => ToBase(strand.SymbolAtHead(head));
        set => strand = strand.WriteSymbolAtHead(head, ToChar(value));
    }

    public Base ReadNextBase()
    {
        var passed = Base;
        if (Direction == ReadingComplexDirection.Forward) head.Move();   // downstream
        else head.MoveReverse();                                         // upstream: the rewind step
        return passed;
    }

    public void AddBase(Base baseToAdd)     => strand = strand.AppendSymbol(ToChar(baseToAdd));
    
    public void Prepend(Base baseToPrepend) => strand = strand.PrependSymbol(ToChar(baseToPrepend));

    // ---- backbone chemistry ----
    public ISequence Cleave()
    {
        var (up, down) = strand.SplitAtHead(head);
        strand = down;
        head.Reset();
        return new FaithfulSequence(up);
    }

    public ISequence Ligate(ISequence other)
    {
        strand = strand + other.ToString();   // join at the tail
        return this;
    }

    // ---- reading-head algorithms (no random access, no length) ----

    // A same-shape backbone of empty sites: walk this strand, laying down one Empty per step.
    public ISequence CreateEmptyBackbone()
    {
        var backbone = new FaithfulSequence("");
        NewRead();
        while (HasNextBase()) { backbone.AddBase(Base.Empty); ReadNextBase(); }
        return backbone;
    }

    // Identity is Watson-Crick closure: lockstep walk, every site must pair, no overhang.
    public bool ClosesCleanlyAgainst(ISequence template)
    {
        NewRead();
        template.NewRead();
        while (HasNextBase() && template.HasNextBase())
        {
            if (!Pairs(Base, template.Base)) return false;
            ReadNextBase();
            template.ReadNextBase();
        }
        return !HasNextBase();   // this strand did not run past the template
    }

    // Does this strand extend beyond the template's end?
    public bool RunsPastEndOf(ISequence template)
    {
        NewRead();
        template.NewRead();
        while (HasNextBase() && template.HasNextBase()) { ReadNextBase(); template.ReadNextBase(); }
        return HasNextBase();
    }

    public bool ContainsNoneOf(Base target)
    {
        NewRead();
        while (HasNextBase())
        {
            if (Base == target) return false;
            ReadNextBase();
        }
        return true;
    } 
 
    static bool Pairs(Base p, Base t)
    {
        if (p == Base.Empty  || t == Base.Empty)  return true;    // abasic pairs anything
        if (p == Base.Border || t == Base.Border) return p == t;  // Border pairs only Border
        return (p == Base.A && t == Base.T) || (p == Base.T && t == Base.A)
            || (p == Base.C && t == Base.G) || (p == Base.G && t == Base.C);
    }

    public override string ToString()
    {
        return strand;
    }
}
