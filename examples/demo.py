"""Run from the project root: python3 -m examples.demo."""

from babai.graph import Graph, graph_isomorphisms
from babai.johnson import johnson_configuration, recognize_johnson


def main():
    print("Babai quasipolynomial graph isomorphism reference implementation")
    cycle = Graph(5, [(i, (i + 1) % 5) for i in range(5)])
    renamed = cycle.relabel((2, 4, 1, 0, 3))
    result = graph_isomorphisms(cycle, renamed)
    assert result is not None
    assert cycle.validates_isomorphism(renamed, result.representative)
    print("An isomorphism of the 5-cycle:", result.representative)
    print("Number of isomorphisms of the 5-cycle:", result.order)

    cycle6 = Graph(6, [(i, (i + 1) % 6) for i in range(6)])
    triangles = Graph(6, [(0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (3, 5)])
    assert graph_isomorphisms(cycle6, triangles) is None
    print("6-cycle versus two triangles: not isomorphic")

    model = recognize_johnson(johnson_configuration(7, 3))
    assert model is not None
    print(f"Johnson scheme: {len(model.subsets)} vertices; recovered {model.atoms} underlying points")


if __name__ == "__main__":
    main()
