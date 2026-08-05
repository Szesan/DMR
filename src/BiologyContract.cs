// The physical interface: what a strand and a medium can do. Constraints in README §2.
using System;
using System.Linq;
/// <summary>Represents a single-stranded nucleic acid polymer processed sequentially by a molecular reading head.</summary>
public interface ISequence 
{
   /// <summary>Travel direction of the reading head, relative to the sense strand.</summary>
   ReadingComplexDirection Direction { get; set; }

   /// <summary>Resets the molecular reading head to the initiation site at the beginning of the strand.</summary>
    void NewRead();

    /// <summary>Checks if the reading head is currently positioned over a valid nucleotide base.</summary>
    bool HasNextBase();

    /// <summary>The base at the reading head's active site.</summary>
    Base Base { get; set; }

    /// <summary>Reads the base at the active site, then advances the head one position downstream.</summary>
    Base ReadNextBase();

    /// <summary>Elongates the polymer by polymerizing and appending a new nucleotide base to the downstream end of the strand.</summary>
    void AddBase(Base baseToAdd);

    /// <summary>Prepends a base at the upstream origin, shifting the strand's register.</summary>
    void Prepend(Base baseToPrepend);

    /// <summary>Synthesizes an empty backbone of this strand's length, all sites abasic.</summary>
    ISequence CreateEmptyBackbone();

    /// <summary>Simulates endonuclease cleavage at the current position of the reading head. This instance is truncated.</summary>
    ISequence Cleave();

    /// <summary>Ligates another strand onto the downstream terminus of this one</summary>
    ISequence Ligate(ISequence other);

    /// <summary>Watson-Crick closure: does this strand base-pair cleanly against the template at every site?</summary>
    bool ClosesCleanlyAgainst(ISequence template);

    /// <summary>Does this strand still have sequence left when the template is spent?</summary>
    bool RunsPastEndOf(ISequence template);

    /// <summary>Does this strand contain NO occurrence of ?</summary>
    bool ContainsNoneOf(Base target);

}

/// <summary>Travel direction of a reading head; Forward is the default, Reverse rewinds upstream.</summary>
public enum ReadingComplexDirection
{
    Forward,
    Reverse,
}

/// <summary>Specifies the chemical and epigenetic identity of a single nucleotide position.</summary>
public enum Base
{
    /// <summary>An abasic (AP) site: the backbone is present but the base is missing.</summary>
    Empty,

    /// <summary>Adenine nucleobase.</summary>
    A,

    /// <summary>Cytosine nucleobase.</summary>
    C,

    /// <summary>Guanine nucleobase.</summary>
    G,

    /// <summary>Thymine nucleobase.</summary>
    T,

    /// <summary>5-hydroxymethyluracil (5hmu). Acts as an epigenetic boundary marker.</summary>
    Border
}

/// <summary>DNA storage medium. It recovers stretches of sequence by their epigenetic boundary motifs</summary>
public interface IGenome
{
    /// <summary>The medium's only retrieval, address-free, based on pattern matching</summary>
    Pool<ISequence> FetchSequencesByFlankMotifs(string leftFlank, string rightFlank);

    /// <summary>As FetchSequencesByFlankMotifs but without the payload cap, so whole clustered copies return intact.</summary>
    Pool<ISequence> FetchClustered(string leftFlank, string rightFlank);

    /// <summary>Clears every stretch bracketed by the given boundary motifs</summary>
    void EliminateByMotif(string leftFlank, string rightFlank);

    /// <summary>Frames each sequence in the given boundary motifs and disperses the copies across distant loci</summary>
    void StoreDistributed(Pool<ISequence> sequences, string leftFlank, string rightFlank);
    
    /// <summary>As but lays the framed copies down together in one contiguous locus</summary>
    void StoreClustered(Pool<ISequence> sequences, string leftFlank, string rightFlank);

    /// <summary>Substitutes the base at a single site but epigenetic boundary marks are immune</summary>
    void ApplyMutation(int position, Base mutatedBase);

    /// <summary>Initiates a new, bare strand — empty and ready to be elongated base by base.</summary>
    ISequence CreateNewSequence();
    /// <summary>Initiates a new, bare double-stranded molecule, both strands empty.</summary>
    DoubleStrand CreateNewDnaSequence();
    /// <summary>Initiates a new RNA transcript — empty and ready to be transcribed.</summary>
    RNA CreateNewRnaSequence();

    void NewRead();

    /// <summary>Runs an opaque job over many independent units at scale, the medium's massively-parallel capacity</summary>
    Pool<TOut> ProcessAtScale<TIn, TOut>(Func<TIn, TOut> job, Pool<TIn> units);

}

