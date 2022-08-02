import itertools
import json
import operator
import time
from collections import deque
from decimal import *  # pylint: disable=wildcard-import
from functools import cached_property, wraps
from math import comb, pi, sqrt
from typing import Callable, Optional

from matplotlib import pyplot as plt
from tqdm.auto import tqdm

from decimal_math import cos, sin

TOLERANCE = Decimal(1e-10)
PRECISION = 10
HIGH_PRECISION = 25

ZERO = Decimal(0)
END = (Decimal(-1), ZERO)
START = (Decimal(1), ZERO)


def fn_timer(fn: Callable) -> Callable:
    """
    Decorates a function to display running time.

    Args:
        fn (function): _description_

    Returns:
        Callable: Decorated function
    """

    @wraps(fn)
    def function_timer(*args, **kwargs):
        start = time.time()
        result = fn(*args, **kwargs)
        end = time.time()
        print(f"Total time running {fn.__name__}: {(end - start):.2f} seconds.")
        return result

    return function_timer


def zeroise(val: Decimal) -> Decimal:
    """
    Returns zero if the absolute value is smaller a given tolerance.
    Otherwise returns the value itself.

    Args:
        val (Optional[float]): _description_

    Returns:
        Optional[float]: _description_
    """
    return val if abs(val) > TOLERANCE else ZERO


class Point:
    """
    Point class that contains the coordinates of a given point.
    """

    def __init__(self, x: Decimal, y: Decimal, precision=PRECISION):
        self.x: Decimal = round(zeroise(x), precision)
        self.y: Decimal = round(zeroise(y), precision)

    @classmethod
    def from_intersection(
        cls, line_a: "LineSegment", line_b: "LineSegment"
    ) -> Optional["Point"]:
        """
        Alternate constructor for Point via intersection of two LineSegments.
        """
        a1, b1, c1 = (
            (line_a.p2.x - line_a.p1.x),
            (line_b.p1.x - line_b.p2.x),
            (line_b.p1.x - line_a.p1.x),
        )
        a2, b2, c2 = (
            (line_a.p2.y - line_a.p1.y),
            (line_b.p1.y - line_b.p2.y),
            (line_b.p1.y - line_a.p1.y),
        )

        d = a1 * b2 - b1 * a2
        d_x = c1 * b2 - b1 * c2
        d_y = a1 * c2 - c1 * a2

        if d != 0:
            s = d_x / d
            t = d_y / d

            if not (0 <= s <= 1 and 0 <= t <= 1):
                return None
        else:
            return None

        pt = cls(
            (1 - s) * line_a.p1.x + s * line_a.p2.x,
            (1 - s) * line_a.p1.y + s * line_a.p2.y,
        )
        return pt.reduced_precision()

    def reduced_precision(self, precision=PRECISION) -> "Point":
        x = round(zeroise(self.x), precision)
        y = round(zeroise(self.y), precision)
        return Point(x, y)

    def value(self, polygon_solver: "PolygonSolver") -> int:
        return polygon_solver.dp[self]

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

    def __init__(self, p1: Point, p2: Point):
        self.p1, self.p2 = sorted([p1, p2], key=operator.attrgetter("x", "y"))

    @property
    def l2(self) -> float:
        return sqrt((self.p1.x - self.p2.x) ** 2 + (self.p1.y - self.p2.y) ** 2)

    @property
    def mid(self) -> Point:
        return Point((self.p1.x + self.p2.x) / 2, (self.p1.y + self.p2.y) / 2)

    def direction(self, pt: Point = None) -> bool:
        """
        Provides the direction of the edge in the graph, by comparing the L2 distances
        (p1-to-end, p2-to-end).

        Caution: This method must only be used when there are no intersection points
            between self.p1 and self.p2. Otherwise, it does not return a useful directional
            label.

        Args:
            pt (Point, optional): True if the direction is pointing towards <pt>.
                                  <pt> should be one of the ends of the line segment.
                                  Defaults to None.

        Returns:
            bool: True if p1-to-p2; False if p2-to-p1.
        """
        direction = (
            LineSegment(self.p1, Point(*END)).l2 > LineSegment(self.p2, Point(*END)).l2
        )

        if pt:
            return direction if pt == self.p2 else (not direction)
        return direction

    def gradient(self) -> tuple[Decimal, Decimal]:
        """
        Creates a tuple used in the visualisation of the direction of an edge.

        Returns:
            tuple[Decimal, Decimal]: A tuple with small values that mathematically represent
                the direction of the edge.
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
                Decimal(multiplier * x_coeff) * delta,
                Decimal(multiplier * y_coeff)
                * delta
                * abs((p2.y - p1.y) / (p2.x - p1.x)),
            )
        except (DivisionByZero, InvalidOperation):
            return (ZERO, multiplier * delta)

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
            **kwargs (optional): Passes any other keyword arguments to self.plot_polygon_graph()
                for convenience.
        """
        self.n = n
        self.plot = plot

        self.polygon = self.create_regular_polygon(n)
        self.line_segments, self.adjacency_list, self.indegrees = self.create_graph()

        # Checks that the generated number of interior points are correct
        self.check_intersections(self.adjacency_list)

        # Solve the graph
        self.dp = self.solve()

        if self.plot:
            if n > 30:
                UserWarning("Plotting for values of n > 30 can be extremely slow.")
            self.plot_graph(**kwargs)

    def find_next_coordinate(self, xy: tuple[Decimal, Decimal], t: Decimal) -> Point:
        """
        Used in polygon creation in cls.create_regular_polygon().

        Args:
            xy (tuple[Decimal, Decimal]): Coordinates(origin vector) of the current point
            t (Decimal): Angle of rotation around the origin

        Returns:
            Point: High precision point
        """
        x, y = xy
        return Point(x * cos(t) - y * sin(t), x * sin(t) + y * cos(t), HIGH_PRECISION)

    def create_regular_polygon(self, n: int) -> list[Point]:
        """
        Defines a regular polygon from its vertices.

        Args:
            n (int): Number of sides of the regular polygon

        Returns:
            list[Point]: Contains the Point classes of the polygon vertices
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
        Returns the list of all basic (unbroken) line segments between polygon vertices.

        Returns:
            list[LineSegment]: Contains basic (unbroken) line segments
        """
        return [LineSegment(a, b) for a, b in itertools.combinations(self.polygon, 2)]

    @cached_property
    def line_segments_with_points(self) -> dict[LineSegment, set[Point]]:
        """
        Returns the mapping of unbroken line segments to all points on the respective line segment.

        Returns:
            dict[LineSegment, set[Point]]:
                Keys: All line segments in the graph
                Values: Every point on a given line segment
        """
        result: dict[LineSegment, set[Point]] = {}
        for line in self.lines:
            result[line] = set(
                (line.p1.reduced_precision(), line.p2.reduced_precision())
            )

        # Intersect every two line segments
        for a, b in itertools.combinations(self.lines, 2):
            point = Point.from_intersection(a, b)

            # Line segments may not intersect
            if point:
                result[a].add(point)
                result[b].add(point)

        return result

    def create_graph(
        self,
    ) -> tuple[set[LineSegment], dict[Point, set[Point]], dict[Point, int]]:
        """
        Creates the adjacency list of the graph by:
        1. Creating all the smaller line segments that make up the edges of the graph
        2. Adding the points to the adjacency list
        """
        print(f"Creating graph for n={self.n}...")

        line_segments: set[LineSegment] = set()
        adjacency_list: dict[Point, set[Point]] = {}
        indegrees: dict[Point, int] = {}

        for _, points in tqdm(self.line_segments_with_points.items()):
            # Create base adjacency list with all keys
            for point in points:
                if point not in adjacency_list:
                    adjacency_list[point] = set()
                    indegrees[point] = 0

            # Sort points by x and y coordinates
            sorted_points = sorted(
                [item.reduced_precision() for item in points],
                key=operator.attrgetter("x", "y"),
            )

            for previous, current in zip(sorted_points, sorted_points[1:]):
                # Create line segment for each successive pair of points on the unbroken line segment
                edge = LineSegment(previous, current)
                line_segments.add(edge)

                # Standardise first and second
                first, second = (current, previous)
                if edge.direction():
                    first, second = second, first

                # Add to adjacency list and increment value of indegrees
                adjacency_list[first].add(second)
                indegrees[second] += 1

        # Checks that the generated number of interior points are correct
        self.check_intersections(adjacency_list)
        return line_segments, adjacency_list, indegrees

    def get_topological_ordering(self) -> list[Point]:
        """
        Returns the topological ordering of the graph.

        Returns:
            list[Point]: List of Points in the graph sorted topologically, starting
                with the start point at (1,0) which has zero indegrees.
        """
        print("Getting topological order...")
        queue = deque([Point(*START)])
        order: list[Point] = []
        adjacency_list = self.adjacency_list.copy()
        indegrees = self.indegrees.copy()

        while queue:
            n = queue.popleft()

            for nb in adjacency_list[n]:
                indegrees[nb] -= 1
                if indegrees[nb] == 0:
                    queue.append(nb)

            del adjacency_list[n]
            order.append(n)

        print("Topological ordering complete.")
        return order

    def solve(self) -> dict[Point, int]:
        """
        Solves the polygons problem progressively (via bottom-up dynamic programming) by:
        1. Topologically sorting the nodes in the graph
        2. Calculating the value of nodes by looking up previously calculated values
        """

        dp = {node: 0 for node in self.adjacency_list}
        dp[Point(*START)] = 1
        order = self.get_topological_ordering()

        print("Summing over the graph...")

        for node in tqdm(order):
            for nb in self.adjacency_list[node]:
                dp[nb] += dp[node]
        return dp

    def result(self) -> int:
        return self.dp[Point(*END)]

    def plot_graph(
        self,
        figsize: int = 8,
        show_lines: bool = True,
        show_points: bool = True,
        show_values: bool = True,
        point_size: float = 1,
        head_width: float = 0.03,
    ) -> None:
        x, y = (
            [point.x for point in self.adjacency_list],
            [point.y for point in self.adjacency_list],
        )

        # Configurations
        limit = 1.1
        offset = Decimal(-0.07)

        _, ax = plt.subplots(figsize=(figsize, figsize))
        ax.set_xlim(-limit, limit)
        ax.set_ylim(-limit, limit)
        ax.set_aspect("equal", adjustable="box")
        ax.grid(axis="both", alpha=0.2)

        # Lines
        if show_lines:
            for line in self.line_segments:
                color = next(
                    ax._get_lines.prop_cycler  # pylint: disable=protected-access
                )["color"]
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
        if show_points:
            ax.scatter(x, y, c="black", s=point_size)

        # Values
        if show_values:
            for pt in self.adjacency_list:
                plt.annotate(pt.value(self), (pt.x + offset, pt.y + offset))

    def check_intersections(self, adjacency_list: dict[Point, set[Point]]) -> None:
        """
        Returns the total number of intersections for an even regular n-sided polygon.
        Formula from https://www.math.uwaterloo.ca/~mrubinst/publications/ngon.pdf

        Raises:
            AssertionError: The number of generated points does not equal the expected
                theoretical number of points.
        """
        print("Checking number of intersection points...")

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

        print(f"n = {n}, total number of points = {len(adjacency_list)}")
        if self.n + interior != len(adjacency_list):
            diff = len(adjacency_list) - int(self.n) - int(interior)
            raise AssertionError(
                f"Expected {int(self.n + interior)} points, got {len(adjacency_list)} points. Difference of {diff}."
            )

    def save_result(self, result: int, filename: str) -> None:
        """
        Saves the obtained result to disk.
        Raises an exception if the current results file has a different value.
        """
        # Read file
        try:
            with open(filename, "r", encoding="utf-8") as f:
                results: dict[str, int] = json.load(f)
        except FileNotFoundError:
            print("Existing results file not found.")
            results = {}

        # Check that the entry is the same
        s = str(self.n)
        extracted_result = results.get(s)
        if extracted_result is None:
            results[s] = result
            with open(filename, "w", encoding="utf-8") as f:
                json.dump(results, f)
            print(f"Added entry for n={self.n} to {filename}.")
        elif extracted_result != result:
            raise AssertionError(
                f"The calculated value for n={self.n} does not equal the existing value!"
            )
        print("Calculated result matches existing result.")


@fn_timer
def main(n: int) -> PolygonSolver:
    polygon_solver = PolygonSolver(n=n, plot=False, figsize=20, show_values=False)
    result = polygon_solver.result()
    polygon_solver.save_result(result, "results.json")
    return polygon_solver


if __name__ == "__main__":
    for i in range(100, 122, 2):
        main(i)
