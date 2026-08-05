// Control model: whole-gene duplication, no validation. See README §1.
using System.Collections.Generic;
using System.Linq;

public class NaiveDuplication : IPreservationModel
{
    readonly IGenome genome;

    public NaiveDuplication(IGenome genome)
    {
        this.genome = genome;
    }

    public IGenome Genome => genome;

    public List<ISequence> RunCycle(int time)
    {
        var copies = genome.FetchClustered(Settings.RecipeLeftFlank, Settings.RecipeRightFlank);
        return copies.ToList();
    }
}
