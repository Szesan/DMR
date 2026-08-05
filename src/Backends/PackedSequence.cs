// Packed medium: same semantics as FaithfulSequence, built to run. See README §4.2.
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

public class PackedSequence : ISimSequence
{
    // 2-bit lane: 32 bases per ulong. Masks: 1 bit per base, 64 per ulong.
    private ulong[] _lane;     // ACGT, 2 bits/base
    private ulong[] _empty;    // Empty (abasic) mask, 1 bit/base
    private ulong[] _border;   // Border (5hmU) mask, 1 bit/base
    private int _startOffset;     // physical base index where this sequence begins
    private int _length;          // active length
    private int _readerPosition;  // reading head, relative to logical 0

    private const char Adenine = 'A';
    private const char Cytosine = 'C';
    private const char Thymine = 'T';
    private const char Guanine = 'G';
    private const char Border = 'h';
    private const char Empty = 'e';

    // lane 2-bit code -> char, for the ACGT case
    private static readonly char[] LaneChar = { Adenine, Cytosine, Guanine, Thymine };

    public PackedSequence(int capacityBases = 16)
    {
        int laneUlongs = Math.Max(1, (capacityBases + 31) / 32);
        int maskUlongs = Math.Max(1, (capacityBases + 63) / 64);
        _lane   = new ulong[laneUlongs];
        _empty  = new ulong[maskUlongs];
        _border = new ulong[maskUlongs];
        _startOffset = 0;
        _length = 0;
        _readerPosition = 0;
    }

    public PackedSequence(string bases) : this(bases.Length)
    {
        for (int i = 0; i < bases.Length; i++)
            AddBase(ToBase(bases[i]));
    }

    // --- storage access (masks first, then lane) ---

    private Base GetBaseAtLogical(int index)
    {
        int p = _startOffset + index;
        int mWord = p >> 6;
        int mBit  = p & 63;
        if (((_border[mWord] >> mBit) & 1UL) != 0) return Base.Border;
        if (((_empty[mWord]  >> mBit) & 1UL) != 0) return Base.Empty;
        int lWord = p >> 5;
        int lShift = (p & 31) << 1;
        ulong code = (_lane[lWord] >> lShift) & 0x3UL;
        return (Base)((int)code + (int)Base.A);   // A,C,G,T contiguous after Empty in the enum
    }

    private void SetBaseAtLogical(int index, Base value)
    {
        int p = _startOffset + index;
        int mWord = p >> 6;
        ulong mBitMask = 1UL << (p & 63);
        int lWord = p >> 5;
        int lShift = (p & 31) << 1;
        ulong laneClear = ~(0x3UL << lShift);

        _empty[mWord]  &= ~mBitMask;
        _border[mWord] &= ~mBitMask;

        switch (value)
        {
            case Base.Empty:
                _empty[mWord] |= mBitMask;
                _lane[lWord] &= laneClear;             // void placeholder 00
                break;
            case Base.Border:
                _border[mWord] |= mBitMask;
                _lane[lWord] &= laneClear;             // void placeholder 00
                break;
            default:
                ulong code = (ulong)((int)value - (int)Base.A) & 0x3UL;
                _lane[lWord] = (_lane[lWord] & laneClear) | (code << lShift);
                break;
        }
    }

    private int LaneCapacity => _lane.Length * 32;
    private int MaskCapacity => _empty.Length * 64;

    private void EnsureCapacityForAppend(int count)
    {
        int needed = _startOffset + _length + count;
        if (needed > LaneCapacity)
        {
            int newUlongs = _lane.Length * 2 + (count / 32) + 1;
            Array.Resize(ref _lane, newUlongs);
        }
        if (needed > MaskCapacity)
        {
            int newUlongs = _empty.Length * 2 + (count / 64) + 1;
            Array.Resize(ref _empty, newUlongs);
            Array.Resize(ref _border, newUlongs);
        }
    }

    private void EnsureCapacityForPrepend(int count)
    {
        if (_startOffset >= count) return;

    // Grow at the front in whole-ulong blocks; shift is a multiple of 64 to keep both planes aligned.
        int laneBlock = Math.Max(4, _lane.Length / 2);
        int shiftBases = ((laneBlock * 32 + 63) / 64) * 64;     // round up to a multiple of 64 bases
        int laneUlongs = shiftBases / 32;
        int maskUlongs = shiftBases / 64;

        var newLane   = new ulong[_lane.Length   + laneUlongs];
        var newEmpty  = new ulong[_empty.Length  + maskUlongs];
        var newBorder = new ulong[_border.Length + maskUlongs];
        Array.Copy(_lane,   0, newLane,   laneUlongs, _lane.Length);
        Array.Copy(_empty,  0, newEmpty,  maskUlongs, _empty.Length);
        Array.Copy(_border, 0, newBorder, maskUlongs, _border.Length);
        _lane = newLane; _empty = newEmpty; _border = newBorder;
        _startOffset += shiftBases;
    }

    // --- ISequence Contract ---

    // Packed backend is a Forward-only speed oracle; it does not model the rewind capability.
    private ReadingComplexDirection _direction = ReadingComplexDirection.Forward;
    public ReadingComplexDirection Direction
    {
        get => _direction;
        set { if (value != ReadingComplexDirection.Forward) throw new NotSupportedException("PackedSequence supports Forward reads only."); _direction = value; }
    }

    public void NewRead() => _readerPosition = 0;

    public bool HasNextBase() => _readerPosition < _length;
    public int Length => _length;

    public bool RunsPastEndOf(ISequence template) => this._length > LengthOf(template);

    private static int LengthOf(ISequence s)
    {
        while (true)
        {
            if (s is ISimSequence sim) return sim.Length;
            if (s is DoubleStrand ds) { s = ds.SenseStrand; continue; }
            if (s is RNA rna) { s = rna.Sequence; continue; }
            return s.ToString().Length;   // foreign backend fallback
        }
    }

    // Empty and Border answer instantly from the mask planes; ACGT falls back to a base scan.
    public bool ContainsNoneOf(Base target)
    {
        if (target == Base.Empty)  return MaskIsClear(_empty);
        if (target == Base.Border) return MaskIsClear(_border);
        // ACGT fallback: scan bases
        for (int i = 0; i < _length; i++)
            if (GetBaseAtLogical(i) == target) return false;
        return true;
    }

    // True iff no mask bit is set across the active region [_startOffset, _startOffset+_length).
    private bool MaskIsClear(ulong[] mask)
    {
        int p = _startOffset;
        int remaining = _length;
        while (remaining > 0)
        {
            int take = Math.Min(64, remaining);
            ulong slice = ReadMask64(mask, p);
            ulong valid = take == 64 ? ~0UL : ((1UL << take) - 1);
            if ((slice & valid) != 0) return false;
            p += take;
            remaining -= take;
        }
        return true;
    }

    // Extract 64 mask bits starting at base position p (across the ulong boundary).
    private static ulong ReadMask64(ulong[] mask, int p)
    {
        int w = p >> 6;
        int off = p & 63;
        ulong lo = mask[w] >> off;
        if (off == 0) return lo;
        ulong hi = (w + 1 < mask.Length) ? mask[w + 1] : 0UL;
        return lo | (hi << (64 - off));
    }

    public Base Base
    {
        get => GetBaseAtLogical(_readerPosition);
        set => SetBaseAtLogical(_readerPosition, value);
    }

    public Base ReadNextBase()
    {
        Base b = GetBaseAtLogical(_readerPosition);
        _readerPosition++;
        return b;
    }

    public void AddBase(Base baseToAdd)
    {
        EnsureCapacityForAppend(1);
        SetBaseAtLogical(_length, baseToAdd);
        _length++;
    }

    public void Prepend(Base baseToPrepend)
    {
        EnsureCapacityForPrepend(1);
        _startOffset--;
        _length++;
        SetBaseAtLogical(0, baseToPrepend);
    }

    public ISequence Cleave()
    {
        var upstream = new PackedSequence(_readerPosition);
        for (int i = 0; i < _readerPosition; i++)
            upstream.AddBase(GetBaseAtLogical(i));

        _startOffset += _readerPosition;
        _length -= _readerPosition;
        _readerPosition = 0;

        return upstream;
    }

    public ISequence Ligate(ISequence other)
    {
        if (other is PackedSequence packedOther)
        {
            EnsureCapacityForAppend(packedOther._length);
            for (int i = 0; i < packedOther._length; i++)
                AddBase(packedOther.GetBaseAtLogical(i));
        }
        else
        {
            string s = other.ToString();
            EnsureCapacityForAppend(s.Length);
            for (int i = 0; i < s.Length; i++)
                AddBase(ToBase(s[i]));
        }
        return this;
    }

    public ISequence CreateEmptyBackbone()
    {
    // Zero-init lane means A, so emptiness is asserted in the mask rather than implicit.
        var seq = new PackedSequence(this._length);
        for (int i = 0; i < this._length; i++)
            seq.AddBase(Base.Empty);
        return seq;
    }

    // --- closure (the SWAR payoff) ---

    public bool ClosesCleanlyAgainst(ISequence template)
    {
        int probeLen = this._length;

    // Unwrap wrapped templates (RNA, DoubleStrand) so the SWAR kernel still applies.
        var pt = AsPackedSequence(template);
        if (pt != null)
        {
            if (probeLen > pt._length) return false;
            return ClosesPacked(pt, probeLen);
        }

    // Fallback for a genuinely non-packed template: correct, just not vectorized.
        int tmplLen = template.ToString().Length;
        if (probeLen > tmplLen) return false;
        template.NewRead();
        for (int i = 0; i < probeLen; i++)
        {
            if (!PairsBase(GetBaseAtLogical(i), template.Base)) return false;
            template.ReadNextBase();
        }
        return true;
    }

    // Unwrap RNA/DoubleStrand to the underlying PackedSequence, or null if not packed underneath.
    private static PackedSequence AsPackedSequence(ISequence s)
    {
        while (true)
        {
            if (s is PackedSequence ps) return ps;
            if (s is RNA rna) { s = rna.Sequence; continue; }
            if (s is DoubleStrand ds) { s = ds.SenseStrand; continue; }
            return null;
        }
    }

    // per-position pairing for the non-packed fallback
    private static bool PairsBase(Base p, Base tb)
    {
        if (p == Base.Empty || tb == Base.Empty) return true;
        if (p == Base.Border || tb == Base.Border) return p == tb;
        return (p == Base.A && tb == Base.T) || (p == Base.T && tb == Base.A)
            || (p == Base.C && tb == Base.G) || (p == Base.G && tb == Base.C);
    }

    // Chunked closure kernel: 32 bases per iteration via whole-word ops, producing one mismatch bit per base.
    // Pairing rules, in precedence order: Empty pairs anything; Border pairs only Border; else ACGT lanes XOR to 0b11.
    // The two strands have independent start offsets, so aligned lane words and mask slices are extracted from each.
    private bool ClosesPacked(PackedSequence t, int n)
    {
        int i = 0;
        while (i < n)
        {
            int chunk = Math.Min(32, n - i);                 // bases handled this iteration
            ulong validMask = chunk == 32 ? ~0U : ((1UL << chunk) - 1);   // low `chunk` bits valid

            int pa = this._startOffset + i;
            int pb = t._startOffset + i;

            // --- aligned 32-base lane words (64 bits each) ---
            ulong laneA = ReadLane32(this._lane, pa);
            ulong laneB = ReadLane32(t._lane,    pb);
            // complementary iff each 2-bit field XORs to 0b11; a field is BAD iff ~xor has any set bit
            ulong notComp = ~(laneA ^ laneB);                // 00 in a field == that field was 0b11 (good)
            // collapse each 2-bit field to 1 bit: bad if either of its two bits is set
            ulong laneBad = (notComp | (notComp >> 1)) & 0x5555555555555555UL;  // bit at even positions
            // pack the even-position bits down to a contiguous 32-bit-per-base mask
            ulong laneMismatch = Pack2to1(laneBad);

            // --- aligned 32-bit mask slices (1 bit per base) ---
            ulong aEmpty  = ReadMask32(this._empty,  pa);
            ulong bEmpty  = ReadMask32(t._empty,     pb);
            ulong aBorder = ReadMask32(this._border, pa);
            ulong bBorder = ReadMask32(t._border,    pb);

            ulong anyEmpty  = aEmpty | bEmpty;               // abasic on either side -> wildcard
            ulong anyBorder = aBorder | bBorder;
            ulong bothBorder = aBorder & bBorder;
            ulong borderBad = anyBorder & ~bothBorder;       // exactly one side Border -> bad

        // Lane is void at Border positions, so Border must not give a lane verdict; Empty wildcards win.
            ulong mismatch = ((laneMismatch & ~anyBorder) | borderBad) & ~anyEmpty;
            mismatch &= validMask;

            if (mismatch != 0) return false;
            i += chunk;
        }
        return true;
    }

    // Extract 32 lane-codes (64 bits) starting at base position `p`, across the ulong boundary.
    private static ulong ReadLane32(ulong[] lane, int p)
    {
        int w = p >> 5;                  // 32 codes per ulong
        int off = (p & 31) << 1;         // bit offset within the word
        ulong lo = lane[w] >> off;
        if (off == 0) return lo;
        ulong hi = (w + 1 < lane.Length) ? lane[w + 1] : 0UL;
        return lo | (hi << (64 - off));
    }

    // Extract 32 mask bits starting at base position `p`, across the 64-bit ulong boundary.
    private static ulong ReadMask32(ulong[] mask, int p)
    {
        int w = p >> 6;                  // 64 bits per ulong
        int off = p & 63;
        ulong lo = mask[w] >> off;
        if (off == 0) return lo & 0xFFFFFFFFUL;
        ulong hi = (w + 1 < mask.Length) ? mask[w + 1] : 0UL;
        return (lo | (hi << (64 - off))) & 0xFFFFFFFFUL;
    }

    // Pack the bits sitting at even positions (0,2,4,…) of a 64-bit word down into the low 32 bits.
    private static ulong Pack2to1(ulong x)
    {
        x &= 0x5555555555555555UL;
        x = (x | (x >> 1))  & 0x3333333333333333UL;
        x = (x | (x >> 2))  & 0x0F0F0F0F0F0F0F0FUL;
        x = (x | (x >> 4))  & 0x00FF00FF00FF00FFUL;
        x = (x | (x >> 8))  & 0x0000FFFF0000FFFFUL;
        x = (x | (x >> 16)) & 0x00000000FFFFFFFFUL;
        return x;
    }

    // --- ISimSequence Contract ---

    public ISimSequence Clone()
    {
        var clone = new PackedSequence(0);
        clone._lane   = (ulong[])this._lane.Clone();
        clone._empty  = (ulong[])this._empty.Clone();
        clone._border = (ulong[])this._border.Clone();
        clone._startOffset = this._startOffset;
        clone._length = this._length;
        clone._readerPosition = this._readerPosition;
        return clone;
    }

    public ISimSequence Concat(ISimSequence other)
    {
        var result = (PackedSequence)this.Clone();
        result.Ligate(other);
        return result;
    }

    public ISimSequence BasesToNewSequence(List<ISequence> sequenceBases)
    {
        var seq = new PackedSequence();
        foreach (ISequence b in sequenceBases)
            seq.Ligate(b);
        return seq;
    }

    public ISimSequence StringToSequence(string sequence) => new PackedSequence(sequence);

    // --- Overrides & Utilities ---

    private Base ToBase(char c) => c switch
    {
        Adenine => Base.A, Cytosine => Base.C, Guanine => Base.G, Thymine => Base.T,
        Border => Base.Border, _ => Base.Empty,
    };

    private char ToChar(Base b) => b switch
    {
        Base.Border => Border,
        Base.Empty  => Empty,
        _ => LaneChar[(int)b - (int)Base.A],
    };

    public override string ToString()
    {
        if (_length == 0) return string.Empty;
        char[] chars = new char[_length];
        for (int i = 0; i < _length; i++)
            chars[i] = ToChar(GetBaseAtLogical(i));
        return new string(chars);
    }

    public bool Equals(ISequence other) => other != null && ToString() == other.ToString();
    public override bool Equals(object obj) => obj is ISequence s && Equals(s);
    public override int GetHashCode() => ToString().GetHashCode();

    public IEnumerator<ISequence> GetEnumerator()
    {
        for (int i = 0; i < _length; i++)
            yield return new PackedSequence(ToChar(GetBaseAtLogical(i)).ToString());
    }

    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();
}
