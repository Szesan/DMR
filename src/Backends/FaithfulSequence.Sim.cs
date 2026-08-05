using System.Collections;
using System.Collections.Generic;

// FaithfulSequence: the ISimSequence surface, kept out of the biological file. See README §3.1.

public partial class FaithfulSequence
{
    public int Length
    {
        get { int n = 0; NewRead(); while (HasNextBase()) { n++; ReadNextBase(); } return n; }
    }

    public ISimSequence Clone() => new FaithfulSequence(strand);
    public ISimSequence Concat(ISimSequence other) => new FaithfulSequence(strand + other.ToString());
    public ISimSequence StringToSequence(string sequence) => new FaithfulSequence(sequence);

    public ISimSequence BasesToNewSequence(List<ISequence> sequenceBases)
    {
        var built = new FaithfulSequence("");
        var pending = new List<ISequence>(sequenceBases);
        while (pending.Count > 0)   // list plumbing (test convenience), not strand traversal
        {
            var first = pending[0]; pending.RemoveAt(0);
            built.Ligate(first);
        }
        return built;
    }

    public bool Equals(ISequence other) => other != null && strand == other.ToString();
    public override bool Equals(object obj) => obj is ISequence s && Equals(s);
    public override int GetHashCode() => strand.GetHashCode();

    public IEnumerator<ISequence> GetEnumerator()
    {
        NewRead();
        while (HasNextBase()) { yield return new FaithfulSequence(ToChar(Base).ToString()); ReadNextBase(); }
    }
    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

}
