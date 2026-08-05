// DMR with the archive refresh switched off. See README §1.
public class DmrNoMaintenance : Dmr
{
    public DmrNoMaintenance(IGenome genome)
        : base(genome) { }

    public override void SwapAndAmplify(Pool<ISequence> genes, Pool<ISequence> maps, Pool<ISequence> recipes)
    {
    }
}
