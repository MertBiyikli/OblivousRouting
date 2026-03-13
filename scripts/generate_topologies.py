import argparse
import itertools
import json
import os
from typing import Any, Dict, Iterable, List, Tuple

import networkx as nx
import numpy as np
import os
import networkx as nx
from typing import Any, Dict, Optional

def _fmt(v: Any) -> str:
    # Keep it simple: numbers as-is, everything else string
    if v is None:
        return ""
    if isinstance(v, bool):
        return "1" if v else "0"
    return str(v)

def write_lgf_arcs(
        G: nx.Graph,
        out_path: str,
        cost_attr: str = "capacity",
        default_cost: float = 1.0,
        make_bidirectional: bool = True,
        node_label_attr: str = "label",
) -> None:
    """
    Export NetworkX undirected graph to LGF using:
      @nodes
      label
      ...

      @arcs
      label   cost
      u v arc_label cost

    By default, each undirected edge becomes two arcs (u->v and v->u).
    Cost is taken from edge attribute `cost_attr` (default: 'capacity').

    Matches the example format you provided.
    """
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

    nodes = sorted(G.nodes())

    # Create stable arc list
    # For undirected G, iterate normalized edges for determinism
    edges = []
    for u, v, data in G.edges(data=True):
        a, b = (u, v) if u <= v else (v, u)
        edges.append((a, b, data))
    edges.sort(key=lambda t: (t[0], t[1]))

    with open(out_path, "w", encoding="utf-8") as f:
        # Nodes section
        f.write("@nodes\n")
        f.write(f"{node_label_attr}\n")
        for n in nodes:
            f.write(f"{_fmt(n)}\n")

        # Arcs section
        f.write("\n@arcs\n")
        f.write("label\tcost\n")

        arc_id = 0
        for u, v, data in edges:
            cost = data.get(cost_attr, default_cost)

            # u -> v
            f.write(f"{_fmt(u)}\t{_fmt(v)}\t{arc_id}\t{_fmt(cost)}\n")
            arc_id += 1

            if make_bidirectional:
                # v -> u
                f.write(f"{_fmt(v)}\t{_fmt(u)}\t{arc_id}\t{_fmt(cost)}\n")
                arc_id += 1


# -------------------------
# LGF writer
# -------------------------

def _fmt_lgf_value(v: Any) -> str:
    """
    LGF is whitespace-separated; strings with spaces should be quoted.
    We'll quote any non-numeric, and escape quotes/backslashes.
    """
    if v is None:
        return ""
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        # Keep stable float formatting
        return repr(float(v))
    if isinstance(v, bool):
        return "1" if v else "0"

    s = str(v)
    needs_quotes = any(ch.isspace() for ch in s) or any(ch in s for ch in ['"', '\\'])
    if needs_quotes:
        s = s.replace("\\", "\\\\").replace('"', '\\"')
        return f"\"{s}\""
    return s


def write_lgf(
        G: nx.Graph,
        out_path: str,
        name: str = "",
        graph_attributes: Dict[str, Any] | None = None,
        node_label_attr: str = "label",
) -> None:
    """
    Write an undirected NetworkX Graph to LEMON .lgf.

    Output structure:

    @nodes
    label <node_attr_1> <node_attr_2> ...
    0     tor          3           ...
    ...

    @edges
    label u v <edge_attr_1> <edge_attr_2> ...
    e0    0 1 1
    ...

    @attributes
    name "..."
    ...
    """
    if graph_attributes is None:
        graph_attributes = {}

    # Determine node order (stable)
    nodes = sorted(G.nodes())

    # Union of node attribute keys (excluding label)
    node_keys = set()
    for n in nodes:
        node_keys.update(G.nodes[n].keys())
    node_keys.discard(node_label_attr)
    node_keys = sorted(node_keys)

    # Union of edge attribute keys
    edge_keys = set()
    for u, v, data in G.edges(data=True):
        edge_keys.update(data.keys())
    edge_keys.discard("label")
    edge_keys = sorted(edge_keys)

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

    with open(out_path, "w", encoding="utf-8") as f:
        # Nodes section
        f.write("@nodes\n")
        header = [node_label_attr] + node_keys
        f.write("\t".join(header) + "\n")
        for n in nodes:
            row = [n] + [G.nodes[n].get(k, "") for k in node_keys]
            f.write("\t".join(_fmt_lgf_value(x) for x in row) + "\n")

        # Edges section
        f.write("\n@edges\n")
        eheader = ["label", "u", "v"] + edge_keys
        f.write("\t".join(eheader) + "\n")

        # Stable edge order
        # For undirected graphs, normalize (min,max) to avoid duplicates
        edges_norm: List[Tuple[int, int, Dict[str, Any]]] = []
        for u, v, data in G.edges(data=True):
            a, b = (u, v) if u <= v else (v, u)
            edges_norm.append((a, b, data))
        edges_norm.sort(key=lambda t: (t[0], t[1]))

        for idx, (u, v, data) in enumerate(edges_norm):
            row = [u, v] + [data.get(k, "") for k in edge_keys]
            f.write("\t".join(_fmt_lgf_value(x) for x in row) + "\n")

        # Attributes section
        f.write("\n@attributes\n")
        if name:
            f.write(f"name\t{_fmt_lgf_value(name)}\n")
        for k, v in graph_attributes.items():
            f.write(f"{k}\t{_fmt_lgf_value(v)}\n")


# -------------------------
# Original generators (kept as-is, except Grid2D bug fix)
# -------------------------

def generateFatClique(parameters):
    numServerPerToR = parameters['numServerPerToR']
    numLocalToR = parameters['numLocalToR']
    numSubblock = parameters['numSubblock']
    numBlock = parameters['numBlock']

    name = 'FatClique-{0}-{1}-{2}-{3}'.format(numServerPerToR, numLocalToR, numSubblock, numBlock)
    G0 = nx.Graph()

    for bl in range(numBlock):
        for sb in range(numSubblock):
            for tr in range(numLocalToR):
                G0.add_node(tr + numLocalToR*sb + (numLocalToR*numSubblock)*bl,
                            type='tor',
                            numServer=numServerPerToR,
                            tor=tr,
                            subblock=sb,
                            block=bl)

    # Add intra local edges
    for bl in range(numBlock):
        for sb in range(numSubblock):
            nodes = [(numLocalToR*numSubblock)*bl + (numLocalToR)*sb + tr for tr in range(numLocalToR)]
            for (i, j) in itertools.combinations(nodes, 2):
                G0.add_edge(i, j, capacity=1)

    # Add subblock edges
    for bl in range(numBlock):
        for tr in range(numLocalToR):
            nodes = [(numLocalToR*numSubblock)*bl + (numLocalToR)*sb + tr for sb in range(numSubblock)]
            for (i, j) in itertools.combinations(nodes, 2):
                G0.add_edge(i, j, capacity=1)

    # Add block edges
    for sb in range(numSubblock):
        for tr in range(numLocalToR):
            nodes = [(numLocalToR*numSubblock)*bl + (numLocalToR)*sb + tr for bl in range(numBlock)]
            for (i, j) in itertools.combinations(nodes, 2):
                G0.add_edge(i, j, capacity=1)

    return G0, name


def generate2LevelClos(parameters):
    numServer = parameters['numServer']

    G0 = nx.Graph()
    lowers = [0, 1, 2, 3]
    uppers = [4, 5]

    numServer = 3
    for n in lowers:
        G0.add_node(n, type='lower', numServer=numServer)
    for n in uppers:
        G0.add_node(n, type='upper', numServer=0)

    for l in lowers:
        for u in uppers:
            G0.add_edge(l, u, capacity=1)

    name = '2LevelClos-3-4-2'
    return G0, name


def generateToyPaper(parameters):
    G0 = nx.Graph()
    lowers = [0, 1, 2, 3]
    uppers = [4, 5]

    for n in lowers:
        if n < 2:
            G0.add_node(n, type='lower', numServer=2)
        else:
            G0.add_node(n, type='lower', numServer=3)

    for n in uppers:
        G0.add_node(n, type='upper', numServer=0)

    for l in lowers:
        for u in uppers:
            if (l, u) in [(2, 5), (3, 5)]:
                G0.add_edge(l, u, capacity=2)
            else:
                G0.add_edge(l, u, capacity=1)

    name = 'ToyPaper'
    return G0, name


def generateToyPaper2(parameters):
    G0 = nx.Graph()
    lowers = [0, 1, 2, 3]
    uppers = [4, 5]

    for n in lowers:
        if n < 2:
            G0.add_node(n, type='lower', numServer=2)
        else:
            G0.add_node(n, type='lower', numServer=3)

    for n in uppers:
        G0.add_node(n, type='upper', numServer=0)

    for l in lowers:
        for u in uppers:
            if (l, u) in [(2, 4), (3, 5)]:
                G0.add_edge(l, u, capacity=2)
            else:
                G0.add_edge(l, u, capacity=1)

    name = 'ToyPaper2'
    return G0, name


def generateFatTreePartial(parameters):
    switchRadix = parameters['switchRadix']
    numAgg = parameters['numAgg']

    assert numAgg <= switchRadix
    numServer = int(switchRadix/2)
    numToRPerAgg = int(switchRadix/2)
    numSwitchPerAgg = int(switchRadix/2)
    numUplinkPerAggSwitch = int(switchRadix/2)
    numSpi = int(numUplinkPerAggSwitch*numSwitchPerAgg*numAgg/switchRadix)

    totToR = numToRPerAgg*numAgg

    def torId(tr, agg):
        return tr + numToRPerAgg*agg

    def anodeId(sw, agg):
        return sw + numSwitchPerAgg*agg + totToR

    def snodeId(spi):
        return spi + numSwitchPerAgg*numAgg + totToR

    G0 = nx.Graph()

    for agg in range(numAgg):
        for tr in range(numToRPerAgg):
            G0.add_node(torId(tr, agg),
                        numServer=numServer,
                        type='tor',
                        tor=tr,
                        agg=agg)
        for sw in range(numSwitchPerAgg):
            G0.add_node(anodeId(sw, agg),
                        numServer=0,
                        type='switch',
                        sw=sw,
                        agg=agg)

    for spi in range(numSpi):
        G0.add_node(snodeId(spi),
                    type='switch',
                    numServer=0,
                    spi=spi)

    for agg in range(numAgg):
        for sw in range(numSwitchPerAgg):
            for tr in range(numToRPerAgg):
                G0.add_edge(anodeId(sw, agg), torId(tr, agg),
                            capacity=1)

    numrep = int(np.floor(numUplinkPerAggSwitch*numSwitchPerAgg*1.0/numSpi))
    numcommonlink = numrep * numSpi
    residue = numUplinkPerAggSwitch*numSwitchPerAgg % numSpi
    for agg in range(numAgg):
        cnt = 0
        for sw in range(numSwitchPerAgg):
            for i in range(numUplinkPerAggSwitch):
                cnt += 1
                if cnt <= numcommonlink:
                    spi = (i + sw*numUplinkPerAggSwitch) % numSpi
                else:
                    spi = (i + sw*numUplinkPerAggSwitch + residue*agg) % numSpi
                aid, sid = anodeId(sw, agg), snodeId(spi)
                if not G0.has_edge(aid, sid):
                    G0.add_edge(aid, sid,
                                capacity=1)
                else:
                    G0[aid][sid]['capacity'] += 1

    for spi in range(numSpi):
        sid = snodeId(spi)
        assert sum(G0[sid][h]['capacity'] for h in G0.neighbors(sid)) == switchRadix

    name = 'FatTree-{0}-{1}'.format(switchRadix, numAgg)
    return G0, name


def generateRing(parameters):
    numServerPerToR = parameters['numServerPerToR']
    numToR = parameters['numToR']
    linkCapacity = parameters['linkCapacity']

    name = 'Ring-{0}-{1}-{2}'.format(numServerPerToR, numToR, linkCapacity)
    G0 = nx.Graph()

    for tr in range(numToR):
        G0.add_node(tr,
                    type='tor',
                    numServer=numServerPerToR,
                    tor=tr)

    for tr in range(numToR):
        if not G0.has_edge(tr, (tr+1) % numToR):
            G0.add_edge(tr, (tr+1) % numToR, capacity=linkCapacity)
        else:
            G0[tr][(tr+1) % numToR]['capacity'] += linkCapacity

    return G0, name


def generateGrid2D(parameters):
    numServerPerToR = parameters['numServerPerToR']
    numRow = parameters['numRow']
    numCol = parameters['numCol']
    linkCapacity = parameters['linkCapacity']  # FIX: was missing in your snippet

    name = 'Grid2D-{0}-{1}-{2}-{3}'.format(numServerPerToR, numRow, numCol, linkCapacity)
    G0 = nx.Graph()

    def nodeid(r, c):
        return c + r*numCol

    for r in range(numRow):
        for c in range(numCol):
            G0.add_node(nodeid(r, c), numServer=numServerPerToR, row=r, col=c)

    for r in range(numRow):
        for c in range(numCol-1):
            G0.add_edge(nodeid(r, c), nodeid(r, c+1), capacity=linkCapacity)

    for c in range(numCol):
        for r in range(numRow-1):
            G0.add_edge(nodeid(r, c), nodeid(r+1, c), capacity=linkCapacity)

    return G0, name


def generateTorus2D(parameters):
    numServerPerToR = parameters['numServerPerToR']
    numRow = parameters['numRow']
    numCol = parameters['numCol']
    linkCapacity = parameters['linkCapacity']

    name = 'Torus2D-{0}-{1}-{2}-{3}'.format(numServerPerToR, numRow, numCol, linkCapacity)
    G0 = nx.Graph()

    def nodeid(r, c):
        return c + r*numCol

    for r in range(numRow):
        for c in range(numCol):
            G0.add_node(nodeid(r, c), numServer=numServerPerToR, row=r, col=c)

    for r in range(numRow):
        for c in range(numCol):
            if not G0.has_edge(nodeid(r, c), nodeid(r, (c+1) % numCol)):
                G0.add_edge(nodeid(r, c), nodeid(r, (c+1) % numCol), capacity=linkCapacity)
            else:
                G0[nodeid(r, c)][nodeid(r, (c+1) % numCol)]['capacity'] += linkCapacity

            if not G0.has_edge(nodeid(r, c), nodeid((r+1) % numRow, c)):
                G0.add_edge(nodeid(r, c), nodeid((r+1) % numRow, c), capacity=linkCapacity)
            else:
                G0[nodeid(r, c)][nodeid((r+1) % numRow, c)]['capacity'] += linkCapacity

    return G0, name


def generateGrid3D(parameters):
    numServerPerToR = parameters['numServerPerToR']
    numRow = parameters['numRow']
    numCol = parameters['numCol']
    numLev = parameters['numLev']
    linkCapacity = parameters['linkCapacity']

    name = 'Grid3D-{0}-{1}-{2}-{3}-{4}'.format(numServerPerToR, numRow, numCol, numLev, linkCapacity)
    G0 = nx.Graph()

    def nodeid(r, c, l):
        return c + r*numCol + l*numRow*numCol

    for r in range(numRow):
        for c in range(numCol):
            for l in range(numLev):
                G0.add_node(nodeid(r, c, l), numServer=numServerPerToR, row=r, col=c, lev=l)

    for r in range(numRow):
        for c in range(numCol):
            for l in range(numLev):
                if r+1 < numRow:
                    G0.add_edge(nodeid(r, c, l), nodeid(r+1, c, l), capacity=linkCapacity)
                if c+1 < numCol:
                    G0.add_edge(nodeid(r, c, l), nodeid(r, c+1, l), capacity=linkCapacity)
                if l+1 < numLev:
                    G0.add_edge(nodeid(r, c, l), nodeid(r, c, l+1), capacity=linkCapacity)

    return G0, name


def generateTorus3D(parameters):
    numServerPerToR = parameters['numServerPerToR']
    numRow = parameters['numRow']
    numCol = parameters['numCol']
    numLev = parameters['numLev']
    linkCapacity = parameters['linkCapacity']

    name = 'Torus3D-{0}-{1}-{2}-{3}-{4}'.format(numServerPerToR, numRow, numCol, numLev, linkCapacity)
    G0 = nx.Graph()

    def nodeid(r, c, l):
        return c + r*numCol + l*numRow*numCol

    for r in range(numRow):
        for c in range(numCol):
            for l in range(numLev):
                G0.add_node(nodeid(r, c, l), numServer=numServerPerToR, row=r, col=c, lev=l)

    for r in range(numRow):
        for c in range(numCol):
            for l in range(numLev):
                if not G0.has_edge(nodeid(r, c, l), nodeid((r+1) % numRow, c, l)):
                    G0.add_edge(nodeid(r, c, l), nodeid((r+1) % numRow, c, l), capacity=linkCapacity)
                else:
                    G0[nodeid(r, c, l)][nodeid((r+1) % numRow, c, l)]['capacity'] += linkCapacity

                if not G0.has_edge(nodeid(r, c, l), nodeid(r, (c+1) % numCol, l)):
                    G0.add_edge(nodeid(r, c, l), nodeid(r, (c+1) % numCol, l), capacity=linkCapacity)
                else:
                    G0[nodeid(r, c, l)][nodeid(r, (c+1) % numCol, l)]['capacity'] += linkCapacity

                if not G0.has_edge(nodeid(r, c, l), nodeid(r, c, (l+1) % numLev)):
                    G0.add_edge(nodeid(r, c, l), nodeid(r, c, (l+1) % numLev), capacity=linkCapacity)
                else:
                    G0[nodeid(r, c, l)][nodeid(r, c, (l+1) % numLev)]['capacity'] += linkCapacity

    return G0, name


def generateDRing(parameters):
    numServerPerToR = parameters['numServerPerToR']
    numToRPerBlock = parameters['numToRPerBlock']
    numBlock = parameters['numBlock']
    linkCapacity = parameters['linkCapacity']

    name = 'DRing-{0}-{1}-{2}'.format(numServerPerToR, numToRPerBlock, numBlock, linkCapacity)
    G0 = nx.Graph()

    def nodeid(tor, block):
        return tor + numToRPerBlock*block

    for block in range(numBlock):
        for tor in range(numToRPerBlock):
            inode = nodeid(tor, block)
            G0.add_node(inode,
                        type='tor',
                        numServer=numServerPerToR,
                        tor=tor,
                        block=block)

    for block in range(numBlock):
        for itor in range(numToRPerBlock):
            inode = nodeid(itor, block)
            for jtor in range(numToRPerBlock):
                jnode = nodeid(jtor, (block+1) % numBlock)
                knode = nodeid(jtor, (block+2) % numBlock)
                G0.add_edge(inode, jnode, capacity=linkCapacity)
                G0.add_edge(inode, knode, capacity=linkCapacity)

    return G0, name


def generateAbstractAugmentExpander(parameters):
    groups = parameters['groups']
    linkCapInt = parameters['linkCapInt']
    linkCapEx = parameters['linkCapEx']

    G0 = nx.Graph()
    name = 'AbstractAugExpander-{0}'.format(len(groups))

    def nodeid(tor, gid):
        offset = 0
        for k in range(gid):
            offset += groups[k]['numToR']
        return tor + offset

    for gid in groups.keys():
        numToR = groups[gid]['numToR']
        numServerPerToR = groups[gid]['numServerPerToR']
        name += '_{0}_{1}'.format(numToR, numServerPerToR)
        for tor in range(numToR):
            nid = nodeid(tor, gid)
            G0.add_node(nid, numServer=numServerPerToR, type='tor',
                        tor=tor, group=gid)

    nodepairid = dict()
    cntnode = G0.number_of_nodes()
    for g1, g2 in itertools.combinations(groups.keys(), 2):
        nid = nodepairid[g1, g2] = nodepairid[g2, g1] = cntnode
        G0.add_node(nid, numServer=0, type='switch', pair=(g1, g2))
        cntnode += 1

    name += '-{0}'.format(len(groups))
    for gid in groups.keys():
        linkcap = linkCapInt[gid]
        name += '_{0}'.format(linkcap)
        if linkcap == 0:
            continue
        nodes = [n for n in G0.nodes() if G0.nodes[n].get('numServer', 0) > 0 and G0.nodes[n].get('group') == gid]
        for n, m in itertools.combinations(nodes, 2):
            G0.add_edge(n, m, capacity=linkcap)

    for g1, g2 in linkCapEx.keys():
        linkcap = linkCapEx[g1, g2]
        name += '-{0}_{1}_{2}'.format(g1, g2, linkcap)
        if linkcap == 0:
            continue
        sid = nodepairid[g1, g2]
        g1nodes = [n for n in G0.nodes() if G0.nodes[n].get('numServer', 0) > 0 and G0.nodes[n].get('group') == g1]
        for n in g1nodes:
            G0.add_edge(n, sid, capacity=linkcap)

    return G0, name


def slimflyModel():
    return [5, 7, 11, 17, 19, 25, 29, 35, 43, 47, 55, 79]


def loadSlimFly(parameters):
    netRadix = parameters['netRadix']
    numServerPerToR = parameters['numServerPerToR']

    sf = open('SlimFly/MMS.' + str(netRadix) + '.adj.txt', 'r')
    numSwitch, numLink = sf.readline().split()
    numSwitch = int(numSwitch)
    numLink = int(numLink)

    G = nx.Graph()
    for s in range(numSwitch):
        G.add_node(s, type='tor', numServer=numServerPerToR)

    for s in range(numSwitch):
        nodelist = sf.readline().split()
        for node in nodelist:
            G.add_edge(s, int(node)-1, capacity=1)
    name = 'Slimfly-{0}-{1}'.format(netRadix, numServerPerToR)
    return G, name


def generateBinaryCube(parameters):
    numServerPerToR = parameters['numServerPerToR']
    numDim = parameters['numDim']
    linkCapacity = parameters['linkCapacity']

    name = 'BCube-{0}-{1}-{2}'.format(numServerPerToR, numDim, linkCapacity)
    G0 = nx.Graph()

    numToR = 2**numDim

    def distance(x, y):
        z = x ^ y
        d = 0
        while z != 0:
            d += z & 1
            z = z >> 1
        return d

    for tr in range(numToR):
        G0.add_node(tr,
                    type='tor',
                    numServer=numServerPerToR,
                    tor=tr)

    for xtr in range(numToR):
        for ytr in range(numToR):
            if distance(xtr, ytr) == 1:
                G0.add_edge(xtr, ytr, capacity=linkCapacity)

    return G0, name


def generateClique(parameters):
    numServerPerToR = parameters['numServerPerToR']
    numToR = parameters['numToR']
    linkCapacity = parameters['linkCapacity']

    name = 'Clique-{0}-{1}-{2}'.format(numServerPerToR, numToR, linkCapacity)
    G0 = nx.Graph()

    tors = list(range(numToR))

    for tr in tors:
        G0.add_node(tr,
                    type='tor',
                    numServer=numServerPerToR,
                    tor=tr)

    for (xtr, ytr) in itertools.combinations(tors, 2):
        G0.add_edge(xtr, ytr, capacity=linkCapacity)

    return G0, name


def generateDRing_80_64(parameters):
    linkCapacity = parameters['linkCapacity']

    name = 'DRing_80_64'
    f = open("dring_80_64.edgelist", "r")
    lines = f.readlines()
    G0 = nx.Graph()
    for inode in range(0, 80):
        if inode in range(0, 52):
            G0.add_node(inode, type='tor', numServer=37)
        else:
            G0.add_node(inode, type='tor', numServer=38)
    for line in lines:
        line = line.replace("\n", "")
        line = line.split("->")
        edge = (int(line[0]), int(line[1]))
        G0.add_edge(edge[0], edge[1], capacity=1)

    return G0, name

def generateRandomRegularExpander(parameters):
    numNodes = parameters['numNodes']
    degree = parameters['degree']
    linkCapacity = parameters['linkCapacity']
    seed = parameters.get('seed', None)

    if numNodes * degree % 2 != 0:
        raise ValueError("numNodes * degree must be even for a d-regular graph.")

    G0 = nx.random_regular_graph(degree, numNodes, seed=seed)

    for u, v in G0.edges():
        G0[u][v]['capacity'] = linkCapacity

    for n in G0.nodes():
        G0.nodes[n]['type'] = 'tor'
        G0.nodes[n]['numServer'] = 0
        G0.nodes[n]['tor'] = n

    name = f'RandomRegularExpander-{numNodes}-{degree}-{linkCapacity}'
    return G0, name

# -------------------------
# CLI glue
# -------------------------

TOPOLOGIES = {
    "fatclique": generateFatClique,
    "2levelclos": generate2LevelClos,
    "toypaper": generateToyPaper,
    "toypaper2": generateToyPaper2,
    "fattree_partial": generateFatTreePartial,
    "ring": generateRing,
    "grid2d": generateGrid2D,
    "torus2d": generateTorus2D,
    "grid3d": generateGrid3D,
    "torus3d": generateTorus3D,
    "dring": generateDRing,
    "abstract_aug_expander": generateAbstractAugmentExpander,
    "slimfly": loadSlimFly,
    "expander": generateRandomRegularExpander,
    "binarycube": generateBinaryCube,
    "clique": generateClique,
    "dring_80_64": generateDRing_80_64,
}

def main():
    ap = argparse.ArgumentParser(description="Generate topology graphs and export to LGF (@arcs format).")

    ap.add_argument("--topology", "-t", choices=sorted(TOPOLOGIES.keys()), required=True)

    # Keep old JSON interface (optional now)
    ap.add_argument("--params", "-p", default=None,
                    help="JSON string or path to JSON file with parameters for the topology.")

    # Output
    ap.add_argument("--out", "-o", required=True, help="Output .lgf file path.")

    # ---- NEW: convenience flags for grids/torus ----
    ap.add_argument("--rows", type=int, default=None, help="(grid/torus) numRow")
    ap.add_argument("--cols", type=int, default=None, help="(grid/torus) numCol")
    ap.add_argument("--lev",  type=int, default=None, help="(3D grid/torus) numLev")
    ap.add_argument("--cap",  type=float, default=None, help="(grid/torus) linkCapacity")
    ap.add_argument("--servers", type=int, default=None, help="(grid/torus) numServerPerToR")

    # LGF options
    ap.add_argument("--directed", action="store_true",
                    help="Write single arc per edge (u->v) instead of bidirectional arcs.")
    ap.add_argument("--cost-attr", default="capacity",
                    help="Edge attribute to export as 'cost' (default: capacity).")

    args = ap.parse_args()

    # Helper: load params JSON if provided
    def load_params_json(p: str) -> dict:
        if os.path.exists(p):
            with open(p, "r", encoding="utf-8") as f:
                return json.load(f)
        return json.loads(p)

    topo = args.topology

    # Decide params source:
    # If user uses convenience flags for grid/torus, build params from flags.
    uses_grid_flags = any(x is not None for x in [args.rows, args.cols, args.lev, args.cap, args.servers])
    is_grid_like = topo in {"grid2d", "torus2d", "grid3d", "torus3d"}

    if is_grid_like and uses_grid_flags:
        # defaults if not provided
        params = {}

        # required for all grid/torus
        if args.rows is None or args.cols is None:
            raise SystemExit("For grid/torus with flags you must pass --rows and --cols.")

        params["numRow"] = args.rows
        params["numCol"] = args.cols
        params["linkCapacity"] = 1.0 if args.cap is None else args.cap
        params["numServerPerToR"] = 0 if args.servers is None else args.servers

        # required for 3D variants
        if topo in {"grid3d", "torus3d"}:
            if args.lev is None:
                raise SystemExit("For grid3d/torus3d with flags you must pass --lev.")
            params["numLev"] = args.lev

    else:
        # legacy mode: must have --params
        if args.params is None:
            raise SystemExit("You must pass --params (JSON string or file) unless using grid/torus flags.")
        params = load_params_json(args.params)

    # Generate
    gen = TOPOLOGIES[topo]
    G, name = gen(params)

    # Export in the exact format your solver expects (@arcs, cost column)
    write_lgf_arcs(
        G,
        args.out,
        cost_attr=args.cost_attr,
        default_cost=1.0,
        make_bidirectional=(not args.directed),
        node_label_attr="label",
    )

    print(f"Wrote {args.out} (name={name}, n={G.number_of_nodes()}, m={G.number_of_edges()})")

if __name__ == "__main__":
    main()