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
const precision = 14
const float64EqualityThreshold = 1e-13

var END Point = Point{-1, 0}
var START Point = Point{1, 0}

// Point
type Point struct {
	x float64
	y float64
}

type LineToPoints = map[LineSegment]mapset.Set[Point]
type PointToLines = map[Point]mapset.Set[LineSegment]
type AdjacencyList = map[Point]mapset.Set[Point]

func roundFloat(val float64) float64 {
	ratio := math.Pow(10, float64(precision))
	return math.Round(val*ratio) / ratio
}

// To allow sorting of points
type ByXY []Point

func (a ByXY) Len() int {
	return len(a)
}

func (a ByXY) Swap(i, j int) {
	a[i], a[j] = a[j], a[i]
}

func (a ByXY) Less(i, j int) bool {
	if a[i].x < a[j].x {
		return true
	} else if a[i].x > a[j].x {
		return false
	} else {
		return a[i].y < a[j].y
	}
}

func almostEqual(a, b float64) bool {
	return math.Abs(a-b) <= float64EqualityThreshold
}

func FromIntersection(L1, L2 LineSegment) (Point, error) {
	// Constructor for Point via the intersection of two line segments
	a1, b1, c1 := (L1.p2.x - L1.p1.x), (L2.p1.x - L2.p2.x), (L2.p1.x - L1.p1.x)
	a2, b2, c2 := (L1.p2.y - L1.p1.y), (L2.p1.y - L2.p2.y), (L2.p1.y - L1.p1.y)

	D := a1*b2 - b1*a2
	Dx := c1*b2 - b1*c2
	Dy := a1*c2 - c1*a2

	var pt Point
	if D != 0 {
		s := Dx / D
		t := Dy / D

		if !(0 <= s && s <= 1 && 0 <= t && t <= 1) {
			return pt, errors.New("no intersection found")
		} else {
			x := (1-s)*L1.p1.x + s*L1.p2.x
			y := (1-s)*L1.p1.y + s*L1.p2.y
			if almostEqual(x, 0) {
				x = 0
			}
			if almostEqual(y, 0) {
				y = 0
			}
			return (Point{roundFloat(x), roundFloat(y)}), nil
		}
	} else {
		return pt, errors.New("no intersection found")
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
	return Point{roundFloat(x), roundFloat(y)}
}

func (ps PolygonSolver) CreateRegularPolygon() []Point {
	// Returns the polygon vertices
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
	// For all combinations of polygon vertices, create line segments
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
	// For all combinations of line segments, create intersection points
	// Returns the map of unbroken line segments to points on the line segments
	// TODO: Make concurrent
	lineToPoints := make(LineToPoints)
	lines := ps.LineSegments()

	// Add both known polygon points onto the line
	for _, line := range lines {
		lineToPoints[line] = mapset.NewSet(line.p1, line.p2)
	}

	// Intersect all line segments
	c := combin.Combinations(len(lines), 2)
	for _, v := range c {
		first, second := lines[v[0]], lines[v[1]]
		point, err := FromIntersection(first, second)
		if err == nil {
			lineToPoints[first].Add(point)
			lineToPoints[second].Add(point)
		}
	}
	return lineToPoints
}

func (ps PolygonSolver) GetFullAdjacencyList(lineToPoints LineToPoints) AdjacencyList {
	// Returns the full adjacency list of nodes and their connected nodes (ignoring direction)
	adjacencyList := make(AdjacencyList)
	for _, points := range lineToPoints {
		// Sort all the points on the unbroken line and get each broken line segment
		sortedPoints := points.ToSlice()
		sort.Stable(ByXY(sortedPoints))

		for i := 1; i < len(sortedPoints); i++ {
			if _, ok := adjacencyList[sortedPoints[i]]; ok {
				adjacencyList[sortedPoints[i]].Add(sortedPoints[i-1])
			} else {
				adjacencyList[sortedPoints[i]] = mapset.NewSet(sortedPoints[i-1])
			}

			if _, ok := adjacencyList[sortedPoints[i-1]]; ok {
				adjacencyList[sortedPoints[i-1]].Add(sortedPoints[i])
			} else {
				adjacencyList[sortedPoints[i-1]] = mapset.NewSet(sortedPoints[i])
			}
		}
	}
	return adjacencyList
}

// Create a channel for every goroutine

// # Approach One: Publish to all relevant channels

// Each vertex has its own channel. When the channel is full, the goroutine kicks off.
// At the end of the goroutine, it writes back to all channels that requires its result.

// The goroutine needs to take in:
// 1. Its input channel where it will read from
// 2. The slice of output channels it needs to write to

// The benefit is that we don't have to topo-sort ahead of time, everything just works as-is.
// You can delegate the responsibility of sending the message to the primary goroutine in a pub-sub model.
// The pub-sub model can allow for symmetry to be handled inside the primary goroutine.

func (ps PolygonSolver) GetPointToChannel(adjacencyList AdjacencyList) map[Point]chan uint {
	pointToChannel := make(map[Point]chan uint)
	for p := range adjacencyList {
		pointToChannel[p] = make(chan uint)
	}
	return pointToChannel
}

func main() {
	n := 6
	ps := PolygonSolver{n}
	lineToPoints := ps.GetLineToPoints()
	adjacencyList := ps.GetFullAdjacencyList(lineToPoints)
	fmt.Println(ps.GetPointToChannel(adjacencyList))
}
