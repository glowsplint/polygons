import itertools
import json
import operator
import time
import warnings
from collections import deque
from decimal import *  # pylint: disable=wildcard-import
from functools import cached_property, wraps
from math import comb, pi
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
        fn (Callable): Function to be decorated

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
        val (Decimal): Decimal to be zeroised
    """
    return val if abs(val) > TOLERANCE else ZERO


def is_multiple(x: int, y: int) -> bool:
    return x % y == 0


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
        Constructor for Point instances via intersection of two LineSegments.
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
        """
        Returns a reduced precision Point instance of itself.

        Args:
            precision (_type_, optional): _description_. Defaults to PRECISION.
        """
        x = round(zeroise(self.x), precision)
        y = round(zeroise(self.y), precision)
        return Point(x, y)

    def value(self, polygon_solver: "PolygonSolver") -> int:
        """
        Returns the calculated value of a node.

        Args:
            polygon_solver (PolygonSolver): PolygonSolver instance
        """
        return polygon_solver.dp[self]

    def order(self, polygon_solver: "PolygonSolver") -> int:
        """
        Returns the topological order of a node, bounded by (0, n-1) inclusive.

        Args:
            polygon_solver (PolygonSolver): PolygonSolver instance
        """
        return polygon_solver.point_ordering[self]

    def __eq__(self, other) -> bool:
        if isinstance(other, Point):
            return not zeroise(self.x - other.x) and not zeroise(self.y - other.y)
        raise TypeError(f"Unable to compare Point and {other}.")

    def __hash__(self) -> int:
        return hash((self.x, self.y))

    def __repr__(self) -> str:
        return f"Point({self.x}, {self.y})"

    def string(self) -> str:
        return f"({self.x},{self.y})"


class LineSegment:

    __slots__ = ["p1", "p2"]

    def __init__(self, p1: Point, p2: Point):
        self.p1, self.p2 = sorted([p1, p2], key=operator.attrgetter("x", "y"))

    @property
    def l2(self) -> Decimal:
        return ((self.p1.x - self.p2.x) ** 2 + (self.p1.y - self.p2.y) ** 2).sqrt()

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
        self.edges, self.adjacency_list, self.indegrees = self.create_graph()
        self.order = self.get_topological_ordering(self.adjacency_list, self.indegrees)
        self.point_ordering, self.graph = self.create_simple_graph(
            self.adjacency_list, self.order
        )

        # Checks that the generated number of interior points are correct
        self.check_all(self.edges, self.adjacency_list)

        # Solve the graph
        self.dp = self.solve(self.adjacency_list, self.indegrees, self.order)

        if self.plot:
            if n > 30:
                UserWarning("Plotting for values of n > 30 can be extremely slow.")
            self.plot_graph(**kwargs)

    def find_next_coordinate(self, xy: tuple[Decimal, Decimal], t: Decimal) -> Point:
        """
        Returns the next clockwise high-precision point from an existing polygon point.

        Args:
            xy (tuple[Decimal, Decimal]): Coordinates (origin vector) of the current point
            t (Decimal): Angle of rotation around the origin

        Returns:
            Point: High precision point
        """
        x, y = xy
        return Point(x * cos(t) - y * sin(t), x * sin(t) + y * cos(t), HIGH_PRECISION)

    def create_regular_polygon(self, n: int) -> list[Point]:
        """
        Returns the vertices of a regular polygon.

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

        edges: set[LineSegment] = set()
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
                edges.add(edge)

                # Standardise first and second
                first, second = (current, previous)
                if edge.direction():
                    first, second = second, first

                # Add to adjacency list and increment value of indegrees
                adjacency_list[first].add(second)
                indegrees[second] += 1

        return edges, adjacency_list, indegrees

    def create_simple_graph(
        self, adjacency_list: dict[Point, set[Point]], topological_order: list[Point]
    ) -> tuple[dict[Point, int], dict[int, set[int]]]:
        """
        Returns the simplified adjacency list that labels each node with their topological order
        from 0 to n-1, 0 being the start point and n-1 being the end point.

        Args:
            adjacency_list (dict[Point, set[Point]]): _description_

        Returns:
            dict[int, set[int]]: _description_
        """
        # Map point to their respective index in the topological order list
        point_ordering: dict[Point, int] = {}
        for i, point in enumerate(topological_order):
            point_ordering[point] = i

        graph: dict[int, set[int]] = {}
        for k, v in adjacency_list.items():
            graph[point_ordering[k]] = set(point_ordering[item] for item in v)
        return point_ordering, graph

    def get_topological_ordering(
        self, adjacency_list: dict[Point, set[Point]], indegrees: dict[Point, int]
    ) -> list[Point]:
        """
        Returns the topological ordering of the graph. Uses Kahn's algorithm.

        Returns:
            list[Point]: List of Points in the graph sorted topologically, starting
                with the start point at (1,0) which has zero indegrees.
        """
        print("Getting topological order...")
        queue = deque([Point(*START)])
        order: list[Point] = []
        adjacency_list_copy = adjacency_list.copy()
        indegrees_copy = indegrees.copy()

        while queue:
            n = queue.popleft()

            for nb in adjacency_list_copy[n]:
                indegrees_copy[nb] -= 1
                if indegrees_copy[nb] == 0:
                    queue.append(nb)

            del adjacency_list_copy[n]
            order.append(n)

        print("Topological ordering complete.")
        return order

    def solve(
        self,
        adjacency_list: dict[Point, set[Point]],
        indegrees: dict[Point, int],
        topological_order: list[Point],
    ) -> dict[Point, int]:
        """
        Solves the polygons problem via bottom-up dynamic programming by:
        1. Topologically sorting the nodes in the graph
        2. Calculating the value of nodes by looking up previously calculated values
        """

        dp = {node: 0 for node in adjacency_list}
        dp[Point(*START)] = 1

        print("Summing over the graph...")

        for node in tqdm(topological_order):
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
        show_order: bool = False,
        point_size: float = 1,
        head_width: float = 0.03,
    ) -> None:
        x, y = (
            [point.x for point in self.adjacency_list],
            [point.y for point in self.adjacency_list],
        )

        # Configurations
        limit = 1.1
        offset = Decimal(-0.05)

        _, ax = plt.subplots(figsize=(figsize, figsize))
        ax.set_xlim(-limit, limit)
        ax.set_ylim(-limit, limit)
        ax.set_aspect("equal", adjustable="box")
        ax.grid(axis="both", alpha=0.2)

        # Warnings
        if show_order and show_values:
            warnings.warn(
                "Both show_order and show_values arguments are set to True. Plot may not be correctly displayed."
            )

        # Lines
        if show_lines:
            for line in self.edges:
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

        # Order
        if show_order:
            for pt in self.adjacency_list:
                plt.annotate(pt.order(self), (pt.x + offset, pt.y + offset))

    def get_theoretical_nodes(self) -> int:
        """
        Returns the theoretical number of nodes in regular n-gon with all diagonals drawn.
        Reference: https://oeis.org/A007569
        """
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
        return int(interior) + n

    def check_nodes(self, adjacency_list: dict[Point, set[Point]]) -> None:
        """
        Checks that the total number of generated nodes are correct.

        Args:
            adjacency_list (dict[Point, set[Point]]): Adjacency list representation of the graph

        Raises:
            ValueError: If the generated and theoretical values do not match
        """
        print("Checking number of intersection points...")
        nodes = self.get_theoretical_nodes()

        if nodes != len(adjacency_list):
            diff = len(adjacency_list) - nodes
            raise ValueError(
                f"Expected {nodes} nodes, got {len(adjacency_list)} nodes. Difference of {diff}."
            )
        print(f"Total number of points = {len(adjacency_list)}")

    def get_theoretical_regions(self) -> int:
        """
        Returns the theoretical number of regions in regular n-gon with all diagonals drawn.
        Reference: https://oeis.org/A007678
        """
        n = self.n
        d = is_multiple

        regions = (
            (
                ((n**4) - 6 * (n**3) + 23 * (n**2) - 42 * n + 24) / 24
                + (-5 * (n**3) + 42 * (n**2) - 40 * n - 48) / 48 * d(n, 2)
                - 3 * n / 4 * d(n, 4)
                + (-53 * (n**2) + 310 * n) / 12 * d(n, 6)
            )
            + 49 * n / 2 * d(n, 12)
            + 32 * n * d(n, 18)
            + 19 * n * d(n, 24)
            - 36 * n * d(n, 30)
            - 50 * n * d(n, 42)
            - 190 * n * d(n, 60)
            - 78 * n * d(n, 84)
            - 48 * n * d(n, 90)
            - 78 * n * d(n, 120)
            - 48 * n * d(n, 210)
        )
        return int(regions)

    def get_theoretical_edges(self) -> int:
        """
        Returns the theoretical number of line segments in regular n-gon with all diagonals drawn.
        Reference: https://oeis.org/A135565
        """
        return self.get_theoretical_nodes() + self.get_theoretical_regions() - 1

    def check_edges(self, edges: set[LineSegment]) -> None:
        """
        Checks that the total number of generated edges are correct.

        Args:
            edges (set[LineSegment]): Edges of the graph

        Raises:
            ValueError: If the generated and theoretical values do not match
        """
        total = self.get_theoretical_edges()
        if total != len(edges):
            raise ValueError(f"expected {total} edges, got {len(edges)} edges")

    def check_all(
        self, edges: set[LineSegment], adjacency_list: dict[Point, set[Point]]
    ) -> None:
        """
        Runs all runtime checks to ensure validity of the generated graph.

        Args:
            edges (set[LineSegment]): Graph edges
            adjacency_list (dict[Point, set[Point]]): Adjacency list representation of the graph
        """
        self.check_nodes(adjacency_list)
        self.check_edges(edges)

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

    def write_adjacency_list_keys_to_disk(self, filename: str) -> None:
        """
        Saves the adjacency list keys to disk.

        Args:
            filename (str): Name of file to write to
        """
        with open(filename, "w", encoding="utf-8") as f:
            points = [point.string() for point in self.adjacency_list]
            json.dump(points, f)
        print(f"Dumped adjacency list keys to {filename}.")


@fn_timer
def main(n: int) -> PolygonSolver:
    polygon_solver = PolygonSolver(
        n=n, plot=True, figsize=20, show_values=False, show_order=True
    )
    result = polygon_solver.result()
    polygon_solver.save_result(result, "results.json")
    return polygon_solver


if __name__ == "__main__":
    solver = main(6)
    # solver.write_adjacency_list_keys_to_disk("adjacency_list.json")
