"""Exercise the installed CLI (BIOMC2_BINDIR) or the current Automake build."""
import csv
import math
import shutil
import os
import random
import subprocess
import tempfile
from collections import Counter
from itertools import permutations
from pathlib import Path

root = Path(__file__).resolve().parents[1]
bin_dir = Path(os.environ.get("BIOMC2_BINDIR", Path.cwd() / "src")).resolve()
example = root / "example"


def run(name, *args):
    return subprocess.check_output([bin_dir / name, *args], text=True,
                                   stderr=subprocess.DEVNULL, timeout=120)


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


def check_assignment_distances():
    rng = random.Random(223)
    trees = []
    for _ in range(8):
        labels = list(range(6))
        rng.shuffle(labels)

        def balanced(tips):
            if len(tips) == 1:
                return tips[0]
            mid = len(tips) // 2
            return (balanced(tips[:mid]), balanced(tips[mid:]))

        trees.append(balanced(labels))

    def split_set(tree):
        result = set()

        def visit(node):
            if isinstance(node, int):
                return 1 << node
            bits = visit(node[0]) | visit(node[1])
            if 2 <= bits.bit_count() <= 4:
                result.add(min(bits, bits ^ 63))
            return bits

        visit(tree)
        return result

    def newick(tree):
        if isinstance(tree, int):
            return f"s{tree}"
        return f"({newick(tree[0])},{newick(tree[1])})"

    for i, a in enumerate(trees):
        for b in trees[i:]:
            left, right = split_set(a), split_set(b)
            rf = len(left ^ right)
            hdist = min(sum(min((x ^ y).bit_count(), 6 - (x ^ y).bit_count())
                            for x, y in zip(left, order)) for order in permutations(right))
            for name, tree in (("a", a), ("b", b)):
                Path(f"{name}.tre").write_text(f"#NEXUS\nBegin trees;\ntree t = {newick(tree)};\nEnd;\n")
            for metric, expected in (("rf", rf), ("hdist", hdist)):
                histogram = dict(tuple(map(int, row.split())) for row in
                                 run("biomc2.sprdist", "a.tre", "b.tre", metric).splitlines())
                assert histogram[expected] == 1 and sum(histogram.values()) == 1

def sampler(tree, alignment, metric, maxdist=8):
    Path("short.ctrl").write_text(f"treefile = {tree}\nseqfile = {alignment}\n"
        f"metric = {metric}\nmaxdist = {maxdist}\nnsegments = 4\nngen = 100\n"
        "warmup = 1\nburnin = 1\nnsamples = 10\n")
    run("biomc2", "short.ctrl")
    trace = Path("post.dist")
    metric = metric if maxdist else "binary"
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
    checker = bin_dir / "likelihood-check"
    if checker.exists():
        run("likelihood-check", example / "paup.tre")
    distances()
    check_assignment_distances()
    for metric in ("dspr", "rf", "hdist"):
        sampler(example / "paup.tre", example / "paup.nexus", metric)
    sampler(example / "paup.tre", example / "paup.nexus", "hdist", 100)
    sampler(example / "paup.tre", example / "paup.nexus", "dspr", 0)
    n, length = 129, 48
    labels = [f"s{i}" for i in range(n)]
    tree_file(Path("large.tre"), labels)
    rng = random.Random(17)
    sequences = ["".join(rng.choices("ACGT", k=length)) for _ in labels]
    Path("large.nex").write_text(f"#NEXUS\nBegin data;\nDimensions ntax={n} nchar={length};\n"
        "Format datatype=DNA interleave;\nMatrix\n" +
        "\n".join(f"{s} {seq}" for s, seq in zip(labels, sequences)) + "\n;\nEnd;\n")
    sampler("large.tre", "large.nex", "dspr")
    if checker.exists():
        run("likelihood-check", "large.tre")
        tree_file(Path("scaling.tre"), [f"s{i}" for i in range(600)])
        run("likelihood-check", "scaling.tre")

    # A parameter word inside a filename must never select that parameter.
    shutil.copyfile(example / "paup.tre", "rate-distance.tre")
    sampler("rate-distance.tre", example / "paup.nexus", "rf")

    # Header-only or inconsistent traces must produce a diagnostic.
    for content in ["", "Nseg 3 Nsamples 1\n", "Nseg 3 Nsamples 1\n4 8\n"]:
        Path("bad.dist").write_text(content)
        result = subprocess.run([bin_dir / "biomc2.summarise", example / "post.tree", "bad.dist"],
                                stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        assert result.returncode > 0

    Path("invalid.ctrl").write_text(f"treefile = {example / 'paup.tre'}\n"
                                   f"seqfile = {example / 'paup.nexus'}\nannealing = invalid\n")
    result = subprocess.run([bin_dir / "biomc2", "invalid.ctrl"],
                            stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    assert result.returncode > 0 and b"ANNEAL" in result.stderr

    # Posterior draws with zero breakpoints still have a valid one-tree mosaic.
    Path("zero.tree").write_text("#NEXUS\nBegin trees;\ntree a = (s0,s1,(s2,s3));\n"
                                "tree b = (s0,s1,(s2,s3));\nEnd;\n")
    Path("zero.dist").write_text("Nseg 2 Nsamples 2 Metric dspr\n4 8\n"
                                "1 -1 0 0 1 1\n2 -1 0 0 1 1\n")
    run("biomc2.summarise", "zero.tree", "zero.dist")
    check_profile(Path("zero.dist"))
    if shutil.which("Rscript"):
        subprocess.run(["Rscript", run_plot if run_plot.exists() else root / "scripts/biomc2.plot"],
                       check=True, stdout=subprocess.DEVNULL)
