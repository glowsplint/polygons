package main

import (
	"errors"
	"fmt"
	"math"

	mapset "github.com/deckarep/golang-set/v2"
	_ "github.com/shopspring/decimal"
	combin "gonum.org/v1/gonum/stat/combin"
)

// Constants
var END Point = Point{-1, 0}
var START Point = Point{1, 0}

// Point
type Point struct {
	x float64
	y float64
}

const float64EqualityThreshold = 1e-15

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

// func (p *Point) ConnectingEdges()

// LineSegment
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

// PolygonSolver
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

func (ps PolygonSolver) GetLineToPoints() map[LineSegment]mapset.Set[Point] {
	lineToPoints := make(map[LineSegment]mapset.Set[Point])
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

func main() {
	// To avoid static code errors for unused vars
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
	fmt.Println((polygonSolver.GetLineToPoints()))
}

// Create all the polygon nodes
// For each polygon node, create a line segment to every other polygon node, and place all the line segments into a slice
// For each line segment,

// Goroutine definition
func Solve() {
	// needs to receive multiple values
	// need to pass the result of this node to downstream nodes

}

// Create one goroutine. This goroutine will spawn more goroutines and pass its result to it
// Infinite loop and select statement, wait until received values on all channels, take sum, then send
