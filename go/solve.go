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
type PointToLines = map[Point]mapset.Set[LineSegment]
type AdjacencyList = map[Point]mapset.Set[Point]

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
	// Alternate constructor for Point via intersection of two LineSegments
	// Returns a pointer to the Point instance
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
			return (Point{x, y}), nil
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
	// Returns the map of unbroken line segments to points on the line segments
	// TODO: Make concurrent
	lineToPoints := make(LineToPoints)
	lines := ps.LineSegments()
	// First add both known polygon points onto the line
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

func (ps PolygonSolver) GetAdjacencyList(lineToPoints LineToPoints) AdjacencyList {
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

func main() {
	// somePoint := Point{3, 10}
	// someLineSegment := LineSegment{Point{0, 0}, Point{1, 2}}
	// fmt.Println(somePoint)
	// fmt.Println(fromIntersection(someLineSegment, someLineSegment))
	// fmt.Println(someLineSegment.L2())
	// fmt.Println(someLineSegment.Mid())
	// fmt.Println(someLineSegment.Direction(someLineSegment.p1))
	// fmt.Println(someLineSegment.OtherEnd(someLineSegment.p1))
	polygonSolver := PolygonSolver{6}
	lineToPoints := polygonSolver.GetLineToPoints()
	// fmt.Println(lineToPoints)
	fmt.Println((polygonSolver.GetAdjacencyList(lineToPoints)))
}

// Create all the polygon nodes
// For each polygon node, create a line segment to every other polygon node, and place all the line segments into a slice
// For each line segment,

// Goroutine definition
// func Solve() {
// needs to receive multiple values
// need to pass the result of this node to downstream nodes

// }

// Approach One: Publish to all relevant channels
// Each vertex has its own buffered channel. When the buffered channel is full, the goroutine kicks off.
// At the end of the goroutine, it writes to all buffered channels that requires its result.
// The goroutine needs to take in:
// 1. Its input buffered channel where it will read from
// 2. The slice of output buffered channels it needs to write to
// The benefit is that we don't have to topo-sort ahead of time, everything just works as-is.
// You can delegate the responsibility of sending the message to the primary goroutine in a pub-sub model.
// The pub-sub model can allow for symmetry to be handled inside the primary goroutine.

// Approach Two:
// 1. Topological sort on the nodes
// 2. Create a queue with these nodes
// 3. Create a pool of worker nodes
// 4. When a task is done, the goroutine sends its output on a single channel.
// 5. The queue processor maintains a single channel and stores the output. It spins up new goroutines when they can be calculated.
// 6. An indegrees array is maintained to trigger the spinning up of new goroutines.

// Optimisations
// 1. We can discard the values of earlier nodes once they are no longer required

// Improvements to be made:
// 1. Using goroutines to speed up processing by saturating available cores
// 2. Using bottom-up dynamic programming and avoiding the recursion limit
// 3. Using symmetry to reduce duplicate calculations

//
