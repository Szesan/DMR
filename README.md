# Distributed Modular Redundancy

A simulation of biological information preservation under continuous mutational damage,
modelled at the granularity of individual base-pair interactions.

The thesis it tests: **an encoding plus a maintenance cycle can preserve genomic information
indefinitely under mutational pressure, without natural selection and without any pristine reference copy.** Everything
lives in the same mutable medium and takes the same damage. Nothing is protected, and no
external template is ever consulted.

---

## 1. Installing and running

### Requirements

.NET SDK 8.0 or later. No external packages, no NuGet dependencies.

### Build and run

```bash
cd src
dotnet run -c Release -- [model] [seed] [ticks] [floor] [--resume]
```

Use `-c Release`. Debug builds run 3–5× slower, because the model reaches the medium
through interface calls that the JIT only devirtualises and inlines under optimisation (§5).
On long runs that is the difference between minutes and hours.

| argument | default | meaning |
|---|---|---|
| model | `fulldmr` | `fulldmr`, `dmr`, or `naive` |
| seed | `100` | seeds genome construction *and* the damage stream |
| ticks | `1000` | damage events to apply |
| floor | `0.5` | stop early if recovered/total falls below this fraction |
| `--resume` | off | continue from a previous run's genome dump |

### The three models

| model | mechanism | outcome |
|---|---|---|
| `fulldmr` | encoding + consensus validation + periodic swap-and-amplify | holds indefinitely |
| `dmr` | encoding + consensus validation, archive never refreshed | long plateau, then collapse |
| `naive` | whole-gene duplication, no validation | degrades from the first ticks |

`dmr` is the informative middle case: it isolates what the encoding alone buys, with the
maintenance cycle switched off. It is the `MaintenanceInterval = ∞` limit of `fulldmr`.

### Output

Written to `runs/`, named `{Model}_seed{N}_maint{N}_report{N}`:

- **`.csv`** — one row per report tick: `tick,recovered,spurious,totalGenes,cycleMs`
- **`.ancestral.txt`** — the pristine genes and their exon bodies, dumped *before any damage*,
  so a reader can independently verify what the run recovers
- **`.genome.txt`** — applied-tick count plus the **state of the medium**. Contrast it against
  `.ancestral.txt` to see the damaged substrate the model was reconstructing from.

  It is rewritten after every maintenance event (`CheckpointOnMaintenance`, default on) and
  again when the run ends, so a long run killed mid-flight stays resumable from its last
  refresh instead of being lost. The write goes to a temp file and is moved into place, so an
  interrupted write cannot truncate the existing checkpoint.

  `--resume` reads this file. Note that resuming replays the damage RNG from tick zero to
  rebuild the identical stream, so resume cost grows with checkpoint age — fast, but not free.

`recovered` counts ancestral genes present in the presented set. `spurious` counts presented
genes that are not ancestral. `cycleMs` times the reconstruction cycle that produced the row —
a diagnostic, not a claim (§9.3).

---

## 2. Constraints

Nothing in the biological path indexes, counts, compares strings, or consults a reference. A
strand is a polymer that a reading head walks along, one base at a time. Every operation
available to the model is something a molecular machine could physically perform.

The commitment is deliberate and expensive, and it is what makes the result meaningful.

### 2.1 The reading head constraint

All strand access goes through a reading head with three operations: read the base at the
current site, step to the next site, and ask whether a site exists.

**Why:** polymerases, ribosomes, and helicases are processive — they bind and track along the
substrate. A molecular machine has a position, not an address space. Random access is the one
assumption that would quietly make the whole problem easy, and the one thing no enzyme can do.

### 2.2 No random access

`IGenome` has exactly
one retrieval operation:

```csharp
Pool<ISequence> FetchSequencesByFlankMotifs(string leftFlank, string rightFlank);
```

Address-free and pattern-matched. You cannot ask for "the sequence at locus 4,812." You
describe the boundary motifs bracketing what you want, and the medium returns every stretch
that matches, wherever it lies. Storage mirrors this: `StoreDistributed` scatters copies across
distant loci, and nothing chooses where.

### 2.3 No loops in the biological sense

The model contains `while` loops, but every one is a head advancing along a strand or a pool
draining one element at a time. There is no iteration over an indexed collection.

### 2.4 No lengths, no counters

`ISequence` has no `Length` member. A head discovers the end of a strand by running off it
(`HasNextBase()` returns false), exactly as a polymerase does.

Counters are excluded likewise. Where a quantity is needed — "make twelve copies" — it is
expressed as a population filling until full, not an integer incrementing to a bound:

The pool *is* the physical population, and asking whether it is full is asking about a
concentration, which a cell can sense in contrast to an abstract counter no molecule performs.

### 2.5 No string identity — only Watson-Crick pairing

Two strands are never compared by string equality in the biological path. The only identity
test available is physical: synthesise the complement of one and attempt to close it against
the other.

```csharp
bool ClosesCleanlyAgainst(ISequence template);
```

Clean closure means every position paired; a single mismatch means failure. This is the *only*
comparison primitive in the model. Consensus validation, type matching in `Balance`, and
seating a probe against a template all reduce to attempted base pairing.

### 2.6 No pristine reference

There is no golden copy anywhere. Every stored sequence sits in the same medium under the same
damage stream. The ancestral genes are written to `.ancestral.txt` at construction for *the
reader's* verification and are never consulted by the model.

This is the constraint that makes the result non-trivial. Any preservation scheme is easy given
an incorruptible reference; the question is what remains possible without one. The answer
demonstrated here: consensus across redundant copies suffices, because health is shared and
decay is idiosyncratic.

---

## 3. Levels of abstraction

Not all code is held to the constraints above. There are three tiers as follows:

**Tier 1 — the model.** `Dmr.cs`, `DmrNoMaintenance.cs`, `NaiveDuplication.cs`, `Enzyme.cs`.
Every constraint in §2 applies without exception. This is the code that should read as biology.

**Tier 2 — the medium.** `Backends/`. Implements what Tier 1 calls. Free to be dense, packed,
and parallel. It implements `ISequence`
and `IGenome` honestly, which the two-backend arrangement (§4) exists to establish.

**Tier 3 — the harness.** `Simulation.cs`, `GenomeFactory.cs`, `Settings.cs`, `Extensions.cs`.
Damage scheduling, scoring, CSV output, genome construction. Makes no biological claim.
`GenomeFactory` is specifically a *construction edge* — it builds the starting genome by fiat,
which is exactly where assertions and preconditions belong so they stay out of the biological
path.

### 3.1 The constraint is structural

The model does not merely decline to use lengths and string equality — it **cannot**. Tier 1
holds only `ISequence`, `DoubleStrand`, `RNA`, and `IGenome` references. Conveniences such as
`Length`, `Concat`, and enumeration live on `ISimSequence`, which the biological path never
touches and cannot reach without a cast. `BiologyContract.cs` contains zero occurrences of the
word `Length`.

---

## 4. The delegate model and the two backends

`ISequence` and `IGenome` are interfaces, and the model calls only those. Two implementations
exist with different roles. 

### 4.1 FaithfulBackend — the demonstration

`FaithfulSequence` honours the reading-head constraint *in its own implementation*. It has a
head position and a direction; it does not index. Cleaving is a split at the head, advancing is
a single step. It is slow, and that is the point: it is the existence proof that every operation
the model performs is achievable under strict processive access.

Its `ISimSequence` members — `Length`, `Clone`, `Concat`, equality, enumeration — are split into
`FaithfulSequence.Sim.cs`, so the biological file contains nothing contradicting the constraint
it demonstrates.

### 4.2 PackedBackend — the performance representation

`PackedSequence` is the same semantics in a form that runs at genome scale. It is permitted to
be dense and to exploit parallelism through `ProcessAtScale`, which fans reconstruction jobs
across cores.

### 4.3 Reasoning for the distinction

The faithful backend proves the operations are physically possible. The packed backend makes
100,000-tick runs feasible. Because both satisfy the same interface and the model holds only the
interface, results are representation-independent — the model cannot tell which one it is
running on.

With only the packed backend, a critic could suspect the mechanism depends on some capability
the packed form quietly provides. With only the faithful backend, the experiments would be too
slow to run.

---

## 5. Release builds 

The interface-delegate pattern is the worst case for unoptimised builds. Every strand operation
is a virtual call through `ISequence` to a small method. Under `-c Release` the JIT observes
that the concrete type is almost always `PackedSequence`, devirtualises, and inlines through it.
Under Debug it does none of that, and every base read remains a real dispatch to a non-inlined
method.

---

## 6. Epigenetic motifs vs. genetic information

The medium distinguishes two kinds of content, and this distinction justifies most of
`Settings.cs`.

**Genetic information** is sequence: `A`, `C`, `G`, `T`. It is what the archive stores, what
degrades, and what must be preserved. It is fully mutable.

**Epigenetic marks** are structural annotations on the medium rather than sequence content:

- `Base.Border` — 5-hydroxymethyluracil, an epigenetic boundary marker
- `Base.Empty` — an abasic (AP) site: backbone present, base absent

Every flank motif is built exclusively from these two symbols:

```csharp
ModulonLeftFlank  = "hhhhhhh";   // 7 Borders
RecipeLeftFlank   = "hhhhhhe";   // 6 Borders + 1 Empty
MapLeftFlank      = "hhhhhee";
MapRightFlank     = "hhhheee";
RecipeRightFlank  = "hhheeee";
ModulonRightFlank = "hheeeee";   // 2 Borders + 5 Empties
```

**Reasoning:** `ApplyMutation` returns early on both `Border` and `Empty` sites, so
epigenetic marks are immune to point mutation. Were the flanks made of ordinary bases, damage
would erase the storage framing itself and the experiment would be measuring the wrong thing.

This is not a convenience granted for the simulation's benefit. 5-hmU is documented at exactly
this role in dinoflagellates, enriched at gene-array boundaries. The framing marks *where things
are*; the sequence carries *what they say*. Only the latter is under mutational pressure.

The six motifs are mutually non-prefixing, so a fetch for one layer cannot accidentally match
another. `SlotFreeMotif` and `SlotOccupied` mark storage availability by the same mechanism.

---

## 7. The process

### 7.1 Clustered damage

Damage is not uniform. Each tick, with probability `ClusterInitiationProbability` (0.25), a
lesion event initiates at a random locus. The centre is hit, and every site within
`ScatterRadius` (15) is independently hit with probability `ScatterProbability` (0.50).

```csharp
void DamageEvent(int length)
{
    if (damageRng.NextDouble() >= Settings.ClusterInitiationProbability) return;
    var centre = damageRng.Next(length);
    Hit(length, centre);
    // ... then every offset within ScatterRadius, each at ScatterProbability
}
```

**Clustered damage:** real damage is not Poisson-uniform. Radiation tracks, oxidative bursts, and
replication-fork collapse produce locally correlated lesions. Uniform damage would be kinder to
a redundancy scheme than reality is, because it spreads thinly across many copies instead of
destroying one thoroughly.

This interacts with distributed storage: because copies are scattered, a
cluster landing on one copy of a module leaves the others untouched.

### 7.2 Distributed storage, and the principle that makes it work

The archive holds three layers, each written by `StoreDistributed`, which frames every payload
in its flank pair and scatters the copies:

| layer | payload | copies |
|---|---|---|
| **roots** | `[modulon] h [exon body]` | x per exon type |
| **maps** | `[archive word][modulon]` | x per type |
| **recipes** | y × `[archive word]` | x per gene |

A *root* is one reusable exon body tagged with its type identifier. A *map* binds an archive
word to a modulon. A *recipe* names the ordered exon types composing one gene.

**The principle: health is shared, decay is idiosyncratic.**

Multiple copies of an exon body were identical when stored. Damage hits them independently, at
different positions. The intact copies therefore still agree with one another, while each
damaged copy disagrees in its own particular way. Agreement is evidence of intactness: not
because any copy is privileged, but because there is one way to be correct and many ways to be
wrong.

This is why no pristine reference is needed. The reference is reconstructed from the population
whenever it is required.

### 7.3 Reconstruction

`RunCycle` reconstructs every gene from the archive, from scratch. It runs the full
reconstruction only on report and maintenance ticks — on other ticks damage accumulates
untouched, so each reported result is a fresh reconstruction against a genome that has taken
every lesion since the last one.

**Step 1 — fetch recipes.** Every surviving recipe copy is retrieved. There is no
deduplication by type: if twelve copies of a gene's recipe survive, twelve reconstruction jobs
run. Only one needs to succeed.

**Step 2 — decode recipes to modulons.** Each recipe is a sequence of archive words; the maps
translate word → modulon. `DecodeRecipeToModulonSequence` walks recipe and map pool in lockstep,
matching by attempted closure.

**Step 3 — fetch candidate roots.** Every root in the archive is retrieved as one pool —
all types together, ~12 copies each. Nothing pre-sorts them by type; which body belongs at which
position is decided by whether it seats, not by a lookup.

**Step 4 — the cartesian approximation.** The candidate space is 12⁵ ≈ 250,000 assemblies per
gene. Enumerating it is infeasible, and a cell does not enumerate. Instead `AssembleGeneVariants`
tries shuffled orderings of the candidate pool (up to `EnrichmentAbundant²`)
seating each candidate against the growing template with `FindsSeatingOn`, which slides the
probe along until it closes cleanly or runs past the end.

This *approximates* the cartesian product by sampling it, which is what a population of
molecules does: many attempts in parallel, most failing, one succeeding.

**Step 5 — two-pass consensus.** Surviving variants are laid down as a tandem array and
validated by consensus, twice:

```csharp
Pool<ISequence> singlyValidated = ConsensusOnTandemArray(LayDownTandemArray(variants));
certified = ConsensusOnTandemArray(LayDownTandemArray(singlyValidated));
```

`ConsensusOnTandemArray` opens the duplex, shifts one strand by exactly one unit against its
neighbour, and attempts closure. Copies that agree close cleanly; copies that differ do not.
Consensus is reduced to a single physical step — no comparator, only geometry.

If no ordering certifies, the whole attempt is retried with the next shuffled ordering. The
gene is either reconstructed or it is not; nothing partial is emitted.

**Justification:** for a mutation to certify it would have to appear at the same position independently in multiple copies, 
twice in a row. Idiosyncratic damage cannot do that. 
The failure mode is "cannot reconstruct," never "reconstructs wrong."

Contrast `naive`, where a damaged copy is indistinguishable from a good one because nothing
arbitrates, and `spurious` climbs monotonically.

The same two-pass consensus runs over the map and recipe layers before any of this, so the
ruleset used for reconstruction is itself validated rather than taken on trust.

### 7.4 Exon validation against the gene template

Consensus in §7.3 certifies *genes*. It says nothing about which archive roots are still intact —
a gene can be reconstructed successfully from an ordering that happened to draw good copies,
while damaged copies of the same exon types sit untouched in the medium.

`ValidateExonsFromGene` closes that gap. It runs before any refresh, and it propagates
gene-level consensus back down to the exon layer.

A certified gene is synthesised into a duplex and its complementary strand becomes a template.
An empty backbone is created alongside it. Then every root in the archive is stripped of its
`[modulon]` tag, leaving a bare exon body, and each body is offered to the template:

```csharp
while ((match = Enzyme.FindsSeatingOn(exon, template)) && !exon.ClosesCleanlyAgainst(senseStrand))
    Enzyme.Shift(exon);

if (match)
{
    Enzyme.Bond(senseStrand, exon);
    validated.Add(exon);
}
```

A body slides along the template until it closes cleanly at a position not already occupied. If
it seats, it is bonded into the growing sense strand and kept. If it never seats, it is dropped
without ceremony — no flag, no report, it simply does not enter the validated pool.

The loop ends when the backbone contains no empty sites: the gene has been fully tiled by intact
bodies. Damaged bodies cannot seat, because seating *is* base pairing against a validated
template.

This is what makes the refresh in §7.5 safe. The exon bodies it copies forward are not merely
*present* in the archive — each has been physically demonstrated to fit inside a gene that
consensus already certified. Validation flows genes → exons, and the refreshed archive is built
only from bodies that passed.

### 7.5 Swap-and-amplify

Reconstruction alone does not correct the archive: degraded roots remain in the medium and are
retrieved again next cycle. Every `MaintenanceInterval` ticks, `SwapAndAmplify` refreshes it.

For each surviving root:

1. Copy the root and cleave off its `[modulon]` tag, keeping both parts.
2. Find a body in the validated pool (§7.4) that pairs with this root's body — this is what
   replaces the possibly-degraded archived body with a demonstrably intact one.
3. Advance the tag under the modulon stencil — the **version bump**.
4. Retire the old lineage: `EliminateByMotif` clears every root carrying the *old* tag.
5. Ligate the advanced tag to the validated body and store the result.
6. Update the maps, so archive words now resolve to the advanced modulons.
7. `Balance` each layer back to full enrichment.

**Why versioning is necessary.** In this model the whole exon pool is refreshed in one
operation. A cell cannot do that — refresh happens in phases, and old and new copies coexist in
the medium during the transition. Without a version marker the compilation machinery could not
tell which generation a root belongs to, and would mix them.

This is what the **modulon / archive word distinction** is for:

- The **modulon** is the *versioned* type tag carried on the root. It advances each generation.
- The **archive word** is the *stable* type name used in recipes. It never changes.
- The **map** binds word → current modulon, and is updated at each refresh.

Recipes therefore never need rewriting: they name types by permanent names, and the map layer
absorbs the versioning. Retiring a whole lineage becomes a single motif elimination against the
old tag.

**The stencil** implements a version bump a molecule could perform:

```csharp
ModulonStencil = "hhheeehhh";   // shoulder | centre | shoulder
```

`AdvanceUnderStencil` advances only the bases exposed through the `Empty` window; the `Border`
shoulders protect the flanks. The tag changes in a bounded, structured way — the type stays
recognisably itself while becoming distinguishable from its predecessor. `MapStencil` does the
same for maps, protecting the 11-base archive-word prefix so only the modulon suffix advances.

---

## 8. Storage budget

The comparison is only meaningful if DMR is not simply given more room. At current settings
(`ExonLength` 90, `GeneCount` 16, 6 exon types, 12× enrichment, 3 naive copies):

| | DMR | naive | ratio |
|---|---|---|---|
| payload bases | 15,384 | 21,600 | **0.71** |
| including flank framing | 20,808 | 22,272 | **0.93** |

DMR is cheaper on both measures — including per-item framing overhead, which counts against it
(288+ separately framed items versus naive's 48). No accounting choice is needed to make the
claim hold.

The medium is 60,000 bases, so both models occupy roughly a third. Note the medium is
*pre-filled with random sequence* and payloads are written into it, so occupancy is not a
meaningful measure of demand.

The cost DMR pays is not storage but computation — elaborate compilation machinery, exactly as
the architecture predicts.

---

## 9. Reading the results

### 9.1 A healthy run

`fulldmr` holds `recovered = totalGenes` and `spurious = 0` at every report tick, indefinitely.
A 100,000-tick run at `MaintenanceInterval` 750 held full recovery across 133 maintenance events
with no deviation in any of 2,000 rows.

### 9.2 The sawtooth

Cycle cost is not flat. It rises as damage accumulates between refreshes and drops sharply at
each maintenance event. The rise is superlinear: as intact witnesses of some exon type thin out,
reconstruction jobs consume more of their 900 shuffled orderings before one certifies.

Peaks are heavy-tailed and driven by *where* damage lands rather than how much has accumulated.
A cluster that removes the last few witnesses of one exon type is expensive; one spread across
well-stocked types is not. Observed peaks have ranged from baseline to 628 s with no loss of
recovery at any point.

### 9.3 What cycle time does and does not mean

Cycle time is a diagnostic, not a failure criterion. The thesis concerns preservation, and a
cycle that takes ten minutes and certifies has preserved the information. A cell has no wall
clock.

It is also ambiguous near collapse: cost is roughly (surviving recipe copies) × (orderings
consumed per job), and damage moves those two terms in opposite directions — fewer surviving
recipes means fewer parallel jobs, which is *cheaper*.

The honest primary statistic is **first-departure tick**: when recovery first leaves
`totalGenes`. It is unambiguous, cheap to observe, and independent of where a floor is set.

### 9.4 The floor drift check

The most informative single measurement for the indefinite-preservation claim is not peak cost
but the *floor* — the cheapest cycle in each window of a long run. If maintenance returned
slightly less redundancy than damage removed, the floor would creep upward across generations.

Across a 100,000-tick run in ten decile windows, the floor sat between 3,643 and 4,547 ms with
no trend. That is the evidence that swap-and-amplify restores the archive fully rather than
partially.

---

## 10. Layout

```
src/
  BiologyContract.cs        ISequence / IGenome — the physical interface, and Base
  SimulationContract.cs     ISimSequence — harness conveniences, invisible to the model
  NucleotideSequence.cs     DoubleStrand and RNA
  Enzyme.cs                 the molecular operations available to a model
  Dmr.cs                    encoding + consensus validation + swap-and-amplify
  DmrNoMaintenance.cs       encoding + validation, no archive refresh
  NaiveDuplication.cs       whole-gene duplication control
  GenomeFactory.cs          builds the starting genome (construction edge)
  Simulation.cs             damage scheduling, scoring, reporting
  Settings.cs               all configuration
  Extensions.cs             Pool<T>, logging, helpers
  Backends/
    FaithfulSequence.cs     strict reading-head medium — the demonstration
    FaithfulSequence.Sim.cs its ISimSequence surface, kept separate
    FaithfulGenome.cs
    PackedSequence.cs       packed medium, same semantics
    PackedGenome.cs
```

---

## 11. Key parameters

| setting | value | role |
|---|---|---|
| `GenomeSize` | 60,000 | medium length in bases |
| `ExonLength` | 90 | one exon body |
| `ExonTypeCount` | 6 | distinct reusable exon types |
| `GeneCount` | 16 | genes composed from those types |
| `ExonsPerGene` | 5 | exons per gene (from 6 types: P(6,5) = 720 orderings) |
| `ModulonLength` | 9 | versioned type tag width |
| `ArchiveWordLength` | 8 | stable type name width |
| `RootEnrichmentRate` | 12 | copies per exon type |
| `IdentifierEnrichmentRate` | 12 | copies per map and per recipe |
| `EnrichmentAbundant` | 30 | squared, gives the 900-ordering retry bound |
| `NaiveCopiesPerGene` | 3 | control model redundancy |
| `MaintenanceInterval` | 600 | ticks between swap-and-amplify |
| `ClusterInitiationProbability` | 0.25 | chance a tick initiates a lesion cluster |
| `ScatterProbability` | 0.50 | chance each site within the radius is hit |
| `ScatterRadius` | 15 | lesion cluster half-width |
| `ShuffleSeed` | 12345 | fixed, so ordering attempts are reproducible |

`MaintenanceInterval` is the parameter deciding whether the cycle outpaces damage.

