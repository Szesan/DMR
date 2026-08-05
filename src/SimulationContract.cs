// Harness-side conveniences. The model never touches these. See README §3.1.
// Testing-convenience layer. 
using System;
using System.Collections.Generic;
public interface ISimSequence : ISequence, IEquatable<ISequence>, IEnumerable<ISequence>
{
    ISimSequence Clone();
    ISimSequence BasesToNewSequence(List<ISequence> bases);
    ISimSequence StringToSequence(string sequence);
    ISimSequence Concat(ISimSequence other);
    int Length { get; }
}

public interface IPreservationModel
{
    /// <summary>Runs one preservation cycle and returns the genes recovered this tick.</summary>
    List<ISequence> RunCycle(int time);

    /// <summary>The storage medium this model preserves into.</summary>
    IGenome Genome { get; }
}