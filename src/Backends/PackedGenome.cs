// Packed 4-bit encoded medium: the performance representation. See README §4.2.


using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

public class PackedGenome : IGenome
{
    private ulong[] _data;
    private int _length;

    private int _distributedCursor;
    private int _clusteredCursor;

    // Fast mappings mirrored exactly
    private const char Adenine = 'A';
    private const char Cytosine = 'C';
    private const char Thymine = 'T';
    private const char Guanine = 'G';
    private const char Border = 'h';
    private const char Empty = 'e';
    
    private static readonly char[] CharMap = { Empty, Adenine, Cytosine, Guanine, Thymine, Border };

    public PackedGenome(ISequence blankMedium)
    {
        string initStr = blankMedium.ToString();
        _length = initStr.Length;
        int ulongs = Math.Max(1, (_length + 15) / 16);
        _data = new ulong[ulongs];
        
        // Fast init
        for (int i = 0; i < _length; i++)
        {
            SetBaseAt(i, ToBase(initStr[i]));
        }

        _distributedCursor = Settings.FirstPosition;
        _clusteredCursor = _length / 2;
    }

    // --- Bitwise Operations (O(1) Access)  --- 

    private Base GetBaseAt(int index)
    {
        if (index >= _length) return Base.Empty;
        return (Base)((_data[index >> 4] >> ((index & 0xF) << 2)) & 0xF);
    }

    private void SetBaseAt(int index, Base value)
    {
        int uIndex = index >> 4;
        int shift = (index & 0xF) << 2;
        _data[uIndex] = (_data[uIndex] & ~(0xFUL << shift)) | ((ulong)value << shift);
    }

    private void EnsureCapacity(int requiredLength)
    {
        int requiredUlongs = (requiredLength + 15) / 16;
        if (requiredUlongs > _data.Length)
        {
            // Amortized doubling identical to List<T> / StringBuilder
            int newSize = Math.Max(requiredUlongs, _data.Length * 2);
            Array.Resize(ref _data, newSize);
        }
    }

    private Base ToBase(char c) => c switch
    {
        Adenine => Base.A, Cytosine => Base.C, Guanine => Base.G, Thymine => Base.T,
        Border => Base.Border, _ => Base.Empty,
    };

    private Base[] ToBaseArray(string str)
    {
        var arr = new Base[str.Length];
        for (int i = 0; i < str.Length; i++) arr[i] = ToBase(str[i]);
        return arr;
    }

    // --- Motif Scanning (String-Free) ---

    private int IndexOf(Base[] pattern, int startIndex)
    {
        if (pattern.Length == 0) return -1;
        int limit = _length - pattern.Length;
        Base first = pattern[0];
        
        for (int i = startIndex; i <= limit; i++)
        {
            if (GetBaseAt(i) == first)
            {
                bool match = true;
                for (int j = 1; j < pattern.Length; j++)
                {
                    if (GetBaseAt(i + j) != pattern[j])
                    {
                        match = false;
                        break;
                    }
                }
                if (match) return i;
            }
        }
        return -1;
    }

    private bool IsLeftFlankAt(int p, Base[] left)
    {
        for (int i = 0; i < left.Length; i++)
        {
            if (GetBaseAt(p + i) != left[i]) return false;
        }
        if (left.Length > 0 && left[0] != Base.Border) return true; // not an h-run
        return !HasBorderAt(p - 1) && !HasBorderAt(p + left.Length);
    }

    private int NextBorderRun(int from, int limit, out int runLength)
    {
        int runStart = -1;
        for (int i = from; i < limit; i++)
        {
            if (GetBaseAt(i) == Base.Border)
            {
                runStart = i;
                break;
            }
        }
        
        if (runStart < 0) 
        {
            runLength = 0;
            return -1;
        }

        int end = runStart;
        while (end < limit && GetBaseAt(end) == Base.Border) end++;
        
        runLength = end - runStart;
        return runStart;
    }

    private bool HasBorderAt(int position) => 
        position >= 0 && position < _length && GetBaseAt(position) == Base.Border;

    // --- IGenome Contract Implementation ---

    public Pool<ISequence> FetchSequencesByFlankMotifs(string leftFlank, string rightFlank)
        => FetchFramed(leftFlank, rightFlank, Settings.MaxPayloadLength);

    // whole-gene retrieval for naive duplication: same scan, no single-module payload cap
    public Pool<ISequence> FetchClustered(string leftFlank, string rightFlank)
        => FetchFramed(leftFlank, rightFlank, int.MaxValue);

    Pool<ISequence> FetchFramed(string leftFlank, string rightFlank, int maxPayload)
    {
        var left = ToBaseArray(leftFlank);
        var right = ToBaseArray(rightFlank);
        var found = new Pool<ISequence>();

        if (left.Length == 0 || right.Length == 0) return found;

        int at = IndexOf(left, 0);
        while (at >= 0)
        {
            int contentStart = at + left.Length;
            int closer = IndexOf(right, contentStart);                     // the next right motif closes the frame
            if (closer >= 0 && closer - contentStart <= maxPayload)
            {
                int payloadLen = closer - contentStart;
                var seq = new PackedSequence(payloadLen);
                for (int i = 0; i < payloadLen; i++) seq.AddBase(GetBaseAt(contentStart + i));
                found.Add(seq);
                at = IndexOf(left, closer + right.Length);
            }
            else at = IndexOf(left, contentStart);
        }
        return found;
    }

    private ISequence ReadFrameUntilExactRun(int contentStart, int rightLength)
    {
        int limit = Math.Min(_length, contentStart + Settings.MaxPayloadLength + rightLength);
        int search = contentStart;
        
        while (search < limit)
        {
            int runStart = NextBorderRun(search, limit, out int runLength);
            if (runStart < 0) break;
            
            if (runLength == rightLength)
            {
                // Directly populate the packed sequence to avoid string allocs
                int payloadLen = runStart - contentStart;
                var seq = new PackedSequence(payloadLen);
                for (int i = 0; i < payloadLen; i++)
                {
                    seq.AddBase(GetBaseAt(contentStart + i));
                }
                return seq;
            }
            search = runStart + runLength;
        }
        return null;
    }

    public void EliminateByMotif(string leftFlank, string rightFlank)
    {
        var left = ToBaseArray(leftFlank);
        var right = ToBaseArray(rightFlank);
        if (left.Length == 0 || right.Length == 0) return;

        // Collect spans first, then overwrite, so editing doesn't disturb the scan.
        var spans = new List<(int start, int length)>();
        int at = IndexOf(left, 0);
        while (at >= 0)
        {
            int contentStart = at + left.Length;
            int closer = IndexOf(right, contentStart);
            if (closer >= 0 && closer - contentStart <= Settings.MaxPayloadLength)
            {
                int end = closer + right.Length;
                spans.Add((at, end - at));
                at = IndexOf(left, end);
            }
            else at = IndexOf(left, contentStart);
        }

        foreach (var (start, length) in spans)                 // overwrite with a neutral base: no h/e survives to re-match
            for (int i = 0; i < length; i++) SetBaseAt(start + i, Base.A);
    }

    private int FrameEndOfExactRun(int contentStart, int rightLength)
    {
        int limit = Math.Min(_length, contentStart + Settings.MaxPayloadLength + rightLength);
        int search = contentStart;
        
        while (search < limit)
        {
            int runStart = NextBorderRun(search, limit, out int runLength);
            if (runStart < 0) break;
            
            if (runLength == rightLength) return runStart + runLength;
            search = runStart + runLength;
        }
        return -1;
    }

    private void WriteAt(int position, string run)
    {
        int required = position + run.Length;
        if (required > _length) EnsureCapacity(required);

        for (int i = 0; i < run.Length; i++)
        {
            SetBaseAt(position + i, ToBase(run[i]));
        }
        if (required > _length) _length = required;
    }

    public void StoreDistributed(Pool<ISequence> sequences, string leftFlank, string rightFlank)
    {
        foreach (var sequence in sequences)
        {
            WriteAt(_distributedCursor, leftFlank + sequence.ToString() + rightFlank);
            _distributedCursor += Settings.ScatterStride;
        }
    }

    public void StoreClustered(Pool<ISequence> sequences, string leftFlank, string rightFlank)
    {
        foreach (var sequence in sequences)
        {
            string run = leftFlank + sequence.ToString() + rightFlank;
            WriteAt(_clusteredCursor, run);
            _clusteredCursor += run.Length;
        }
    }



    public void ApplyMutation(int position, Base mutatedBase)
    {
        if (position < 0 || position >= _length) return;
        if (GetBaseAt(position) == Base.Border) return;
        if (GetBaseAt(position) == Base.Empty) return;
        SetBaseAt(position, mutatedBase);
    }

    public ISequence CreateNewSequence() => new PackedSequence();
    
    public DoubleStrand CreateNewDnaSequence() => 
        new DoubleStrand(new PackedSequence(), new PackedSequence());
    
    public RNA CreateNewRnaSequence() => new RNA(new PackedSequence());

    public Pool<TOut> ProcessAtScale<TIn, TOut>(Func<TIn, TOut> job, Pool<TIn> units)
    {
        return new Pool<TOut>(units.AsParallel().AsOrdered().Select(job));
    }


    public void NewRead()
    {
        _distributedCursor = Settings.FirstPosition;
        _clusteredCursor = _length / 2;
    }

    public override string ToString()
    {
        if (_length == 0) return string.Empty;
        
        // Single block allocation for maximum speed
        char[] chars = new char[_length];
        for (int i = 0; i < _length; i++)
        {
            chars[i] = CharMap[(int)GetBaseAt(i)];
        }
        return new string(chars);
    }
}