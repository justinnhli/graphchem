#!/usr/bin/env python3


from collections import defaultdict


class Timeline:

    # acyclic; cycles of temperature/pressure must have explicit time dependence

    def __init__(self):
        self.temp_pres = []
        self.graph = defaultdict(set)
        self.changed = False
        self.predecessors = defaultdict(set)
        self.topo_order = []

    def __getitem__(self, tpid):
        return self.temp_pres[tpid]

    def add(self, temp_pres):
        self.temp_pres.append(temp_pres)
        return len(self.temp_pres) - 1

    def add_transition(self, tpid1, tpid2):
        self.graph[tpid1].add(tpid2)
        self.changed = True

    def is_predecessor(self, tpid1, tpid2):
        if self.changed:
            self._build_cache()
        return tpid1 in self.predecessors.get(tpid2, set())

    def _build_cache(self):
        srcs = defaultdict(set)
        for src, dsts in self.graph.items():
            for dst in dsts:
                srcs[dst].add(src)
        self.predecessors = defaultdict(set)
        self.topo_order = []
        queue = list(self.graph.keys() - srcs.keys())
        visited = set()
        while queue:
            node = queue.pop()
            visited.add(node)
            self.topo_order.append(node)
            self.predecessors[node].update([node], *(self.predecessors[src] for src in srcs[node]))
            queue.extend(
                successor for successor in self.graph.get(node, set())
                if srcs[successor] <= visited
            )
        self.changed = False

    def convergence_points(self, *tpids):
        if self.changed:
            self._build_cache()
        srcs = set(tpids)
        results = set()
        for node in self.topo_order:
            preds = self.predecessors[node]
            if srcs <= preds and not preds.intersection(results):
                results.add(node)
        return results


def main():
    dsts = {
        'A': {'B', 'C'},
        'B': {'D'},
        'C': {'E', 'G'},
        'D': {'E', 'F'},
        'E': {'F'},
    }
    timeline = Timeline()
    for parent, children in dsts.items():
        timeline.add(parent)
        for child in children:
            timeline.add(child)
    for parent, children in dsts.items():
        parent_id = timeline.temp_pres.index(parent)
        for child in children:
            timeline.add_transition(parent_id, timeline.temp_pres.index(child))

    conv = timeline.convergence_points(
        timeline.temp_pres.index('C'),
        #timeline.temp_pres.index('F'),
    )
    print([timeline.temp_pres[i] for i in conv])

if __name__ == '__main__':
    main()
