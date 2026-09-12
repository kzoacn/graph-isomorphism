import pytest

from babai.blocks import maximal_block_system
from babai.group import PermutationGroup
from babai.permutation import from_cycles
from babai.string_iso import GiantActionRequired, LuksSolver
from babai.tower import BlockTower, refines


def test_tower_restriction_preserves_the_inherited_quotient():
    tower = BlockTower(8).extend(((0, 1), (2, 3), (4, 5), (6, 7))).extend(((0, 1, 2, 3), (4, 5, 6, 7)))
    child = tower.restrict((7, 4, 6, 5))
    assert child.top == ((0, 2), (1, 3))
    assert child.degree == 4
    assert child.restrict((2, 0)).top == ((0,), (1,))
    assert BlockTower(0).restrict(()).top == ()
    with pytest.raises(ValueError, match="coarsenings"):
        tower.extend(((0, 2, 4, 6), (1, 3, 5, 7)))


def test_maximal_blocks_coarsen_the_existing_level():
    cyclic = PermutationGroup(12, [from_cycles(12, range(12))])
    initial = ((0, 4, 8), (1, 5, 9), (2, 6, 10), (3, 7, 11))
    top = maximal_block_system(cyclic, initial)
    assert refines(initial, top)
    assert top == ((0, 2, 4, 6, 8, 10), (1, 3, 5, 7, 9, 11))


def test_solver_carries_tower_into_small_orbit_calls():
    group = PermutationGroup(8, [from_cycles(8, (0, 1)), from_cycles(8, (0, 2, 4, 6), (1, 3, 5, 7)), from_cycles(8, (0, 2), (1, 3))])
    initial = BlockTower(8).extend(((0, 1), (2, 3), (4, 5), (6, 7)))

    class Recorder(LuksSolver):
        def __init__(self):
            super().__init__()
            self.seen = []

        def _solve_instance(self, group, source, target):
            self.seen.append(self.current_tower)
            return super()._solve_instance(group, source, target)

    solver = Recorder()
    result = solver.solve(group, "00112233", "22330011", tower=initial)
    assert result is not None and result.order == 16
    assert solver.seen[0] == initial
    assert any(tower.degree == 2 for tower in solver.seen)
    assert solver.current_tower is None


def test_tower_stack_is_restored_after_a_missing_continuation(monkeypatch):
    monkeypatch.setattr("babai.string_iso.small_quotient_bound", lambda _: 1)
    cyclic = PermutationGroup(5, [from_cycles(5, range(5))])
    solver = LuksSolver()
    with pytest.raises(GiantActionRequired) as exception:
        solver.solve(cyclic, "ababa", "aabab")
    assert exception.value.context.tower.degree == 5
    assert solver.current_tower is None
