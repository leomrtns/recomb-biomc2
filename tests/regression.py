"""Exercise the installed CLI (BIOMC2_BINDIR) or the current Automake build."""
import csv
import math
import shutil
import os
import random
import subprocess
import tempfile
from collections import Counter
from pathlib import Path

root = Path(__file__).resolve().parents[1]
bin_dir = Path(os.environ.get("BIOMC2_BINDIR", Path.cwd() / "src")).resolve()
example = root / "example"


def run(name, *args):
    return subprocess.check_output([bin_dir / name, *args], text=True)


def check_profile(trace):
    lines = trace.read_text().splitlines()
    positions = list(map(int, lines[1].split()))
    samples = [line.split("|") for line in lines[2:]]
    assert all(math.isfinite(float(s[0].split()[1])) for s in samples)
    rows = list(csv.DictReader(Path("breakpoint-profile.tsv").open(), delimiter="\t"))
    for group in ("all", "modal"):
        series = [r for r in rows if r["group"] == group]
        assert {int(r["position"]) for r in series} == set(positions[:-1])
        mode = int(series[0]["modal_count"])
        counts = Counter(len(s) - 1 for s in samples)
        assert counts[mode] == max(counts.values())
        keep = [s for s in samples if group == "all" or len(s) - 1 == mode]
        hits, distances = Counter(), Counter()
        for sample in keep:
            for pair in sample[1:]:
                index, d = map(int, pair.split())
                hits[positions[index]] += 1
                distances[positions[index]] += d
        for r in series:
            pos = int(r["position"])
            assert (int(r["weight"]), int(r["break_weight"]), int(r["distance_sum"])) == (
                len(keep), hits[pos], distances[pos])


def tree_file(path, labels):
    def balanced(tips):
        if len(tips) == 1:
            return tips[0]
        mid = len(tips) // 2
        return f"({balanced(tips[:mid])},{balanced(tips[mid:])})"
    # Quartet with a large common clade. The NNI involves the last two taxon IDs.
    clade = balanced(labels[:-3])
    a, b, c = labels[-3:]
    path.write_text(f"#NEXUS\nBegin trees;\ntree t = [&U] ({clade},{a},({b},{c}));\nEnd;\n")


def distances():
    for n in (4, 63, 64, 65, 129, 257):
        labels = [f"s{i}" for i in range(n)]
        tree_file(Path("a.tre"), labels)
        labels[-2], labels[-3] = labels[-3], labels[-2]
        tree_file(Path("b.tre"), labels)
        for metric, expected in (("rf", 2), ("hdist", 2), ("dspr", 1)):
            for other, d in (("a.tre", 0), ("b.tre", expected)):
                hist = run("biomc2.sprdist", "a.tre", other, metric)
                assert dict(tuple(map(int, row.split())) for row in hist.splitlines())[d] == 1, (metric, n, d)


def sampler(tree, alignment, metric):
    Path("short.ctrl").write_text(f"treefile = {tree}\nseqfile = {alignment}\n"
        f"metric = {metric}\nmaxdist = 8\nnsegments = 4\nngen = 100\n"
        "warmup = 1\nburnin = 1\nnsamples = 10\n")
    run("biomc2", "short.ctrl")
    trace = Path("post.dist")
    assert f"Metric {metric}" in trace.read_text().splitlines()[0]
    run("biomc2.summarise", "post.tree", trace)
    check_profile(trace)


with tempfile.TemporaryDirectory(dir=Path.cwd()) as tmp:
    os.chdir(tmp)
    run("biomc2.summarise", example / "post.tree", example / "post.dist")
    check_profile(example / "post.dist")
    if shutil.which("Rscript"):
        run_plot = bin_dir / "biomc2.plot"
        subprocess.run(["Rscript", run_plot if run_plot.exists() else root / "scripts/biomc2.plot"], check=True)
        assert Path("recomb_freq.pdf").read_bytes().startswith(b"%PDF-")
    run("biomc2.distance", example / "paup.tre", "2")
    assert Path("spr.trees").read_text().rstrip().endswith("end;")
    distances()
    for metric in ("dspr", "rf", "hdist"):
        sampler(example / "paup.tre", example / "paup.nexus", metric)
    n, length = 129, 48
    labels = [f"s{i}" for i in range(n)]
    tree_file(Path("large.tre"), labels)
    rng = random.Random(17)
    sequences = ["".join(rng.choices("ACGT", k=length)) for _ in labels]
    Path("large.nex").write_text(f"#NEXUS\nBegin data;\nDimensions ntax={n} nchar={length};\n"
        "Format datatype=DNA interleave;\nMatrix\n" +
        "\n".join(f"{s} {seq}" for s, seq in zip(labels, sequences)) + "\n;\nEnd;\n")
    sampler("large.tre", "large.nex", "dspr")
