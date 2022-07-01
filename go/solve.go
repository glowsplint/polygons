package main

import (
	"errors"
	"fmt"
	"math"
	"sort"

	mapset "github.com/deckarep/golang-set/v2"
	_ "github.com/shopspring/decimal"
	combin "gonum.org/v1/gonum/stat/combin"
)

// Constants
const float64EqualityThreshold = 1e-15

var END Point = Point{-1, 0}
var START Point = Point{1, 0}

// Point
type Point struct {
	x float64
	y float64
}

type LineToPoints = map[LineSegment]mapset.Set[Point]

type ByXY []Point

func (a ByXY) Len() int {
	return len(a)
}

func (a ByXY) Swap(i, j int) {
	a[i], a[j] = a[j], a[i]
}

func (a ByXY) Less(i, j int) bool {
	return a[i].x < a[j].x && a[i].y < a[j].y
}

func almostEqual(a, b float64) bool {
	return math.Abs(a-b) <= float64EqualityThreshold
}

func FromIntersection(L1, L2 LineSegment) (Point, error) {
	// Alternate constructor for Point via intersection of two LineSegments
	// Returns a pointer to the Point instance
	a1, b1, c1 := (L1.p2.x - L1.p1.x), (L2.p1.x - L2.p2.x), (L2.p1.x - L1.p1.x)
	a2, b2, c2 := (L1.p2.y - L1.p1.y), (L2.p1.y - L2.p2.y), (L2.p1.y - L1.p1.y)

	D := a1*b2 - b1*a2
	Dx := c1*b2 - b1*c2
	Dy := a1*c2 - c1*a2

	if D != 0 {
		s := Dx / D
		t := Dy / D

		if !(0 <= s && s <= 1 && 0 <= t && t <= 1) {
			return Point{0, 0}, errors.New("no intersection found")
		} else {
			return (Point{(1-s)*L1.p1.x + s*L1.p2.x, (1-s)*L1.p1.y + s*L1.p2.y}), nil
		}
	} else {
		return Point{0, 0}, errors.New("no intersection found")
	}
}

type LineSegment struct {
	p1 Point
	p2 Point
}

func (ls *LineSegment) L2() float64 {
	return math.Sqrt(math.Pow(ls.p1.x-ls.p2.x, 2) + math.Pow(ls.p1.y-ls.p2.y, 2))
}

func (ls *LineSegment) Mid() Point {
	return Point{(ls.p1.x + ls.p2.x) / 2, (ls.p1.y + ls.p2.y) / 2}
}

func (ls *LineSegment) Direction(pt Point) bool {
	direction := (&LineSegment{ls.p1, END}).L2() > (&LineSegment{ls.p2, END}).L2()
	if pt == ls.p2 {
		return direction
	} else {
		return !direction
	}
}

func (ls *LineSegment) OtherEnd(pt Point) Point {
	if pt == ls.p1 {
		return ls.p2
	} else {
		return ls.p1
	}
}

type PolygonSolver struct {
	n int
}

func FindNextCoordinate(pt Point, theta float64) Point {
	x := pt.x*math.Cos(theta) - pt.y*math.Sin(theta)
	y := pt.x*math.Sin(theta) + pt.y*math.Cos(theta)
	if almostEqual(x, 0) {
		x = 0
	}
	if almostEqual(y, 0) {
		y = 0
	}
	return Point{x, y}
}

func (ps PolygonSolver) CreateRegularPolygon() []Point {
	theta := 2.0 / float64(ps.n) * math.Pi
	polygon := make([]Point, 0)
	previous := START

	for i := 0; i < ps.n; i++ {
		polygon = append(polygon, previous)
		previous = FindNextCoordinate(previous, theta)
	}
	return polygon
}

func (ps PolygonSolver) LineSegments() []LineSegment {
	c := combin.Combinations(ps.n, 2)
	p := ps.CreateRegularPolygon()
	lineSegments := make([]LineSegment, len(c))
	for i, v := range c {
		start, end := p[v[0]], p[v[1]]
		lineSegments[i] = LineSegment{start, end}
	}
	return lineSegments
}

func (ps PolygonSolver) GetLineToPoints() LineToPoints {
	// TODO: Make concurrent
	lineToPoints := make(LineToPoints)
	lines := ps.LineSegments()
	for _, line := range lines {
		lineToPoints[line] = mapset.NewSet(line.p1, line.p2)
	}

	c := combin.Combinations(len(lines), 2)
	for _, v := range c {
		start, end := lines[v[0]], lines[v[1]]
		point, err := FromIntersection(start, end)

		if err == nil {
			lineToPoints[start].Add(point)
			lineToPoints[end].Add(point)
		}
	}
	return lineToPoints
}

func (ps PolygonSolver) GetAdjacencyList() map[Point]LineSegment {
	// Returns the adjacency list of the graph for every vertex (which has both polygon vertices and
	// intersection points).
	var adjacencyList map[Point]LineSegment
	return adjacencyList
}

func (ps PolygonSolver) GetBrokenLineSegments(lineToPoints LineToPoints) mapset.Set[LineSegment] {
	brokenLineSegments := mapset.NewSet[LineSegment]()
	for line, points := range lineToPoints {
		// A line segment with only two points on it, is already the smallest possible line segment
		if points.Cardinality() == 2 {
			brokenLineSegments.Add(line)
			continue
		}

		// Sort all the points on the unbroken line and get each broken line segment
		sortedPoints := points.ToSlice()
		sort.Sort(ByXY(sortedPoints))

		for i := 1; i < len(sortedPoints); i++ {
			if sortedPoints[i] == sortedPoints[i-1] {
				continue
			}
			edge := LineSegment{sortedPoints[i], sortedPoints[i-1]}
			brokenLineSegments.Add(edge)
		}
	}
	return brokenLineSegments
}

func main() {
	// somePoint := Point{3, 10}
	// someLineSegment := LineSegment{Point{0, 0}, Point{1, 2}}
	// fmt.Println(somePoint)
	// fmt.Println(fromIntersection(someLineSegment, someLineSegment))
	// fmt.Println(someLineSegment.L2())
	// fmt.Println(someLineSegment.Mid())
	// fmt.Println(someLineSegment.Direction(someLineSegment.p1))
	// fmt.Println(someLineSegment.OtherEnd(someLineSegment.p1))
	// required := mapset.NewSet[string]()
	// required.Add("biology")
	// fmt.Println(required)
	polygonSolver := PolygonSolver{6}
	lineToPoints := polygonSolver.GetLineToPoints()
	fmt.Println((polygonSolver.GetBrokenLineSegments(lineToPoints)))
}

// Create all the polygon nodes
// For each polygon node, create a line segment to every other polygon node, and place all the line segments into a slice
// For each line segment,

// Goroutine definition
// func Solve() {
// needs to receive multiple values
// need to pass the result of this node to downstream nodes

// }

// Create one goroutine. This goroutine will spawn more goroutines and pass its result to it
// Infinite loop and select statement, wait until received values on all channels, take sum, then send

// Improvements to be made:
// 1. Using goroutines to speed up processing by saturating available cores
// 2. Using bottom-up dynamic programming to avoid hitting the recursion limit
