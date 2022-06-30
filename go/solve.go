package main

import (
	"errors"
	"fmt"
	"math"

	_ "github.com/shopspring/decimal"
)

// Constants
var END Point = Point{-1, 0}
var START Point = Point{1, 0}

// Point
type Point struct {
	x float64
	y float64
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

	// self.polygon = self.create_regular_polygon(n)
	// self.line_segments = None
	// self.point_ls = None
	// self.create_line_segments()
}

func FindNextCoordinate(pt Point, theta float64) Point {
	return Point{pt.x*math.Cos(theta) - pt.y*math.Sin(theta), pt.x*math.Sin(theta) + pt.y*math.Cos(theta)}
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

// func (ps PolygonSolver) Lines() []LineSegment {
// 	list := combin.Combinations(ps.n, 2)
// 	allLineSegments := make([]LineSegment, 0)
// 	for _, v := range list {
// 		allLineSegments = append(allLineSegments, v)
// 	}
// 	return allLineSegments
// }

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
	polygonSolver := PolygonSolver{4}
	fmt.Println(polygonSolver.CreateRegularPolygon())
}

// Define PolygonSolver
// Define main()

// Create all the polygon nodes
// For each polygon node, create a line segment to every other polygon node, and place all the line segments into a slice
// For each line segment,
