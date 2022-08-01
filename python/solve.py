import itertools
import json
import operator
import sys
import time
from decimal import *
from functools import cached_property, lru_cache, wraps
from math import comb, pi, sqrt

from matplotlib import pyplot as plt
from tqdm.auto import tqdm

from decimal_math import cos, sin

TOLERANCE = Decimal(1e-10)
PRECISION = 13

END = (-1, 0)
START = (1, 0)

sys.setrecursionlimit(30000)


def fn_timer(function):
    @wraps(function)
    def function_timer(*args, **kwargs):
        t0 = time.time()
        result = function(*args, **kwargs)
        t1 = time.time()
        print(f"Total time running {function.__name__}: {(t1 - t0):.2f} seconds.")
        return result

    return function_timer


def zeroise(a: float or None) -> bool:
    if a is None:
        return None
    return a if abs(a) > TOLERANCE else 0


class Point:
    def __init__(self, x: float, y: float, precision=PRECISION):
        self.x = round(zeroise(Decimal(x)), precision)
        self.y = round(zeroise(Decimal(y)), precision)

    @classmethod
    def from_intersection(
        cls, L1: "LineSegment" = None, L2: "LineSegment" = None
    ) -> "Point":
        """
        Alternate constructor for Point via intersection of two LineSegments.
        """
        a1, b1, c1 = (L1.p2.x - L1.p1.x), (L2.p1.x - L2.p2.x), (L2.p1.x - L1.p1.x)
        a2, b2, c2 = (L1.p2.y - L1.p1.y), (L2.p1.y - L2.p2.y), (L2.p1.y - L1.p1.y)

        D = a1 * b2 - b1 * a2
        Dx = c1 * b2 - b1 * c2
        Dy = a1 * c2 - c1 * a2

        if D != 0:
            s = Dx / D
            t = Dy / D

            if not (0 <= s <= 1 and 0 <= t <= 1):
                return None
        else:
            return None

        return cls((1 - s) * L1.p1.x + s * L1.p2.x, (1 - s) * L1.p1.y + s * L1.p2.y)

    def reduced_precision(self, precision=PRECISION) -> "Point":
        x = round(zeroise(Decimal(self.x)), precision)
        y = round(zeroise(Decimal(self.y)), precision)
        return Point(x, y)

    def connecting_edges(self, polygon_solver: "PolygonSolver") -> set["LineSegment"]:
        """
        Provides the adjacent edges to the node in the graph.

        Returns:
            set[LineSegment]: Adjacent edges to the node in the graph.
        """
        return polygon_solver.point_ls[self]

    @lru_cache(maxsize=None)
    def dependent_edges(self, polygon_solver: "PolygonSolver") -> set["LineSegment"]:
        """
        Provides the edges that we depend on to calculate self.value.
        These are the edges that are directed towards the self node.

        Returns:
            set[LineSegment]: Edges that we depend on to calculate self.value
        """
        return set(
            edge
            for edge in self.connecting_edges(polygon_solver)
            if edge.direction(self)
        )

    @lru_cache(maxsize=None)
    def value(self, polygon_solver: "PolygonSolver") -> int:
        """
        Recursively calculates the number of paths we can take from the start node to this node by summing over
            the nodes of dependent edges.

        Args:
            polygon_solver (PolygonSolver): PolygonSolver instance

        Returns:
            int: Number of paths we can take from the start node to this node
        """
        if not self.dependent_edges(polygon_solver):
            return 1
        return sum(
            edge.other_end(self).value(polygon_solver)
            for edge in self.dependent_edges(polygon_solver)
        )

    def show_dependencies(self, polygon_solver, n: int) -> list["Point"]:
        """
        Returns the n-th level of dependent nodes in the graph.
        Primary tool used in debugging RecursionError where nodes have cyclical dependence.

        Returns:
            list[Point]: list of Point nodes at the n-th level of dependence
        """
        current = set([self])
        for _ in range(n):
            current = set(
                edge.other_end(point)
                for point in current
                for edge in point.dependent_edges(polygon_solver)
            )
        return sorted(list(current), key=operator.attrgetter("x", "y"))

    def __eq__(self, other) -> bool:
        if isinstance(other, Point):
            return not zeroise(self.x - other.x) and not zeroise(self.y - other.y)
        raise TypeError(f"Unable to compare Point and {other}.")

    def __hash__(self) -> int:
        return hash((self.x, self.y))

    def __repr__(self) -> str:
        return f"Point({self.x}, {self.y})"


class LineSegment:

    __slots__ = ["p1", "p2"]

    def __init__(self, p1: "Point", p2: "Point"):
        self.p1, self.p2 = sorted([p1, p2], key=operator.attrgetter("x", "y"))

    @property
    def L2(self) -> float:
        return sqrt((self.p1.x - self.p2.x) ** 2 + (self.p1.y - self.p2.y) ** 2)

    @property
    def mid(self) -> "Point":
        return Point((self.p1.x + self.p2.x) / 2, (self.p1.y + self.p2.y) / 2)

    def direction(self, pt: "Point" = None) -> bool:
        """
        Provides the direction of the edge in the graph, by comparing the L2 distances (p1-to-end, p2-to-end).

        Caution: This method must only be used when there are no intersection points between self.p1 and self.p2.
            Otherwise, it does not return a useful directional label.

        Args:
            pt (Point, optional): True if the direction is pointing towards <pt>.
                                  <pt> should be one of the ends of the line segment. Defaults to None.

        Returns:
            bool: True if p1-to-p2; False if p2-to-p1.
        """
        direction = (
            LineSegment(self.p1, Point(*END)).L2 > LineSegment(self.p2, Point(*END)).L2
        )

        if pt:
            return direction if pt == self.p2 else (not direction)
        return direction

    def gradient(self) -> tuple[Decimal, Decimal]:
        """
        Creates a tuple used in the visualisation of the direction of an edge.

        Returns:
            tuple[Decimal, Decimal]: A tuple with small values that mathematically represent the direction of the edge.
        """
        p1, p2 = self.p1, self.p2
        multiplier = 1 if self.direction() else -1
        x_coeff = 0 if p1.x == p2.x else 1
        delta = Decimal(0.01)

        if p2.y > p1.y:
            y_coeff = 1
        elif p1.y == p2.y:
            y_coeff = 0
        else:
            y_coeff = -1

        try:
            return (
                multiplier * x_coeff * delta,
                multiplier * y_coeff * delta * abs((p2.y - p1.y) / (p2.x - p1.x)),
            )
        except (DivisionByZero, InvalidOperation):
            return (0, multiplier * delta)

    def other_end(self, p: "Point") -> "Point":
        """
        Returns:
            Point: Other end of the line segment
        """
        if p == self.p1:
            return self.p2
        elif p == self.p2:
            return self.p1

    def __eq__(self, other) -> bool:
        if isinstance(other, LineSegment):
            return self.p1 == other.p1 and self.p2 == other.p2
        return False

    def __hash__(self) -> int:
        return hash((self.p1, self.p2))

    def __repr__(self) -> str:
        return f"LineSegment({self.p1}, {self.p2})"


class PolygonSolver:
    def __init__(self, n: int, plot: bool = False, **kwargs):
        """
        Generates the entire graph of edges and points upon instantiation.

        Args:
            n (int): Number of sides of the regular polygon
            plot (bool, optional): Plots the graph if true. Defaults to False.
            **kwargs (optional): Passes any other keyword arguments to self.plot_polygon_graph() for convenience.
        """
        self.n = n
        self.polygon = self.create_regular_polygon(n)

        self.line_segments = None
        self.point_ls = None
        self.create_line_segments()

        if plot:
            self.plot_graph(**kwargs)

    @classmethod
    def find_next_coordinate(self, xy: tuple[Decimal, Decimal], t: Decimal) -> "Point":
        """
        Used in polygon creation in cls.create_regular_polygon().

        Args:
            xy (tuple[Decimal, Decimal]): Coordinates(origin vector) of the current point
            t (Decimal): Angle of rotation around the origin

        Returns:
            Point: High precision point
        """
        x, y = xy
        return Point(x * cos(t) - y * sin(t), x * sin(t) + y * cos(t), 25)

    @classmethod
    def create_regular_polygon(self, n: int) -> list[tuple]:
        """
        Defines a regular polygon from its vertices.

        Args:
            n (int): Number of sides of the regular polygon

        Returns:
            list[tuple]: Contains the Point classes of the polygon vertices
        """
        theta = Decimal(2 / n * pi)
        polygon = []
        previous = Point(*START)

        for _ in range(n):
            polygon.append(previous)
            previous = self.find_next_coordinate((previous.x, previous.y), theta)

        return polygon

    @cached_property
    def lines(self) -> list[LineSegment]:
        """
        Returns the list of all basic line segments between polygon vertices.
        """
        return [LineSegment(a, b) for a, b in itertools.combinations(self.polygon, 2)]

    @cached_property
    def line_to_points(self) -> dict[LineSegment, set[Point]]:
        """
        Returns the mapping of unbroken line segments to all points on the respective line segment.

        Returns:
            dict[LineSegment, set[Point]]:
                Keys: All line segments in the graph
                Values: Every point on a given line segment
        """
        line_to_points = dict()
        for line in self.lines:
            line_to_points[line] = set([line.p1, line.p2])

        for a, b in itertools.combinations(self.lines, 2):
            point = Point.from_intersection(a, b)

            if point is not None:
                line_to_points[a].add(point)
                line_to_points[b].add(point)

        return line_to_points

    def create_line_segments(self) -> None:
        """
        Creates all the smaller line segments that make up the edges of the graph.
        """
        line_segments = set()
        point_ls = dict()

        print("Creating line segments...")

        for ls, p in tqdm(self.line_to_points.items()):
            if len(p) == 2:
                line_segments.add(ls)
                continue

            sorted_points = sorted(
                [item.reduced_precision() for item in p],
                key=operator.attrgetter("x", "y"),
            )

            for previous, current in zip(sorted_points, sorted_points[1:]):
                edge = LineSegment(previous, current)
                line_segments.add(edge)
                for item in (previous, current):
                    if item not in point_ls:
                        point_ls[item] = set([edge])
                    else:
                        point_ls[item].add(edge)

        self.line_segments = line_segments
        self.point_ls = point_ls

        # Checks that the generated number of interior points are correct
        self.check_intersections()
        return None

    def plot_graph(
        self,
        figsize: int = 8,
        lines: bool = True,
        points: bool = True,
        values: bool = True,
        point_size: float = 1,
        head_width: float = 0.03,
    ) -> None:
        x, y = (
            [point.x for point in self.point_ls],
            [point.y for point in self.point_ls],
        )

        # Configurations
        limit = 1.1
        offset = -0.07

        _, ax = plt.subplots(figsize=(figsize, figsize))
        ax.set_xlim(-limit, limit)
        ax.set_ylim(-limit, limit)
        ax.set_aspect("equal", adjustable="box")
        ax.grid(axis="both", alpha=0.2)

        # Lines
        if lines:
            for line in self.line_segments:
                color = next(ax._get_lines.prop_cycler)["color"]
                ax.plot([line.p1.x, line.p2.x], [line.p1.y, line.p2.y], color=color)
                gradient = line.gradient()
                plt.arrow(
                    x=line.mid.x,
                    y=line.mid.y,
                    dx=gradient[0],
                    dy=gradient[1],
                    lw=0,
                    color=color,
                    length_includes_head=True,
                    head_width=head_width,
                )

        # Points
        if points:
            ax.scatter(x, y, c="black", s=point_size)

        # Values
        if values:
            for pt in self.point_ls:
                plt.annotate(pt.value(self), (pt.x + offset, pt.y + offset))

    def check_intersections(self) -> None:
        """
        Returns the total number of intersections for an even regular n-sided polygon.
        Formula from https://www.math.uwaterloo.ca/~mrubinst/publications/ngon.pdf

        Raises:
            AssertionError: The number of generated points does not equal the expected
                theoretical number of points.
        """

        def is_multiple(x: int, y: int) -> bool:
            return x % y == 0

        n = self.n
        d = is_multiple
        interior = (
            comb(n, 4)
            + (-5 * n**3 + 45 * n**2 - 70 * n + 24) / 24 * d(n, 2)
            - (3 * n / 2) * d(n, 4)
            + (-45 * n**2 + 262 * n) / 6 * d(n, 6)
            + 42 * n * d(n, 12)
            + 60 * n * d(n, 18)
            + 35 * n * d(n, 24)
            - 38 * n * d(n, 30)
            - 82 * n * d(n, 42)
            - 330 * n * d(n, 60)
            - 144 * n * d(n, 84)
            - 96 * n * d(n, 90)
            - 144 * n * d(n, 120)
            - 96 * n * d(n, 210)
        )

        print(f"n = {n}, total number of points = {len(self.point_ls)}")
        if self.n + interior != len(self.point_ls):
            diff = len(self.point_ls) - int(self.n) - int(interior)
            raise AssertionError(
                f"Expected {int(self.n + interior)} points, got {len(self.point_ls)} points. Difference of {diff}."
            )

    def solve(self) -> int:
        return sorted(self.point_ls.keys(), key=operator.attrgetter("x", "y"))[0].value(
            self
        )

    def progressive_solve(self) -> int:
        """
        Intended as an alternative solver that allows us to progressively solve the graph.
        Allows us to see a progress bar with a marginal performance hit compared to PolygonSolver.solve().
        """
        iterable = sorted(
            self.point_ls.keys(), key=operator.attrgetter("x", "y"), reverse=True
        )
        print("Solving progressively...")
        for item in tqdm(iterable, total=len(iterable)):
            item.value(self)
        return iterable[-1].value(self)

    def save_result(self, value, filename: str) -> None:
        """
        Saves the obtained result to disk.
        Raises an exception if the current results file has a different value.
        """
        # Read file
        with open(filename, "r") as f:
            results: dict[str, int] = json.load(f)

        # Check that the entry is the same
        s = str(self.n)
        extracted_result = results.get(s)
        if extracted_result is None:
            results[s] = value
            with open(filename, "w") as f:
                json.dump(results, f)
            print(f"Added entry for n={self.n} to {filename}.")

        if extracted_result != value:
            raise AssertionError(
                f"The calculated value for n={self.n} does not equal the existing value!"
            )
        else:
            print("Calculated result matches existing result.")


@fn_timer
def main():
    solver = PolygonSolver(n=6, plot=True, figsize=20, values=False)
    final = solver.progressive_solve()
    print(final)
    solver.save_result(final, "results.json")
    return solver, final


if __name__ == "__main__":
    solver, final = main()
