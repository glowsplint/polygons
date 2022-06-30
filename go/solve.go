package main

import (
	"errors"
	"fmt"
	"math"
	"strings"

	dec "github.com/shopspring/decimal"
)

// Constants
var ZERO dec.Decimal = dec.NewFromInt(0)
var ONE dec.Decimal = dec.NewFromInt(1)
var TWO dec.Decimal = dec.NewFromInt(2)

var END Point = Point{ONE.Neg(), ZERO}
var START Point = Point{ONE, ZERO}

// Point
type Point struct {
	x dec.Decimal
	y dec.Decimal
}

func (p Point) ToString() string {
	return strings.Join([]string{p.x.String(), p.y.String()}, ", ")
}

func FromIntersection(L1, L2 LineSegment) (Point, error) {
	// Alternate constructor for Point via intersection of two LineSegments
	// Returns a pointer to the Point instance
	a1, b1, c1 := (L1.p2.x.Sub(L1.p1.x)), (L2.p1.x.Sub(L2.p2.x)), (L2.p1.x.Sub(L1.p1.x))
	a2, b2, c2 := (L1.p2.y.Sub(L1.p1.y)), (L2.p1.y.Sub(L2.p2.y)), (L2.p1.y.Sub(L1.p1.y))

	D := a1.Mul(b2).Sub(b1.Mul(a2))
	Dx := c1.Mul(b2).Sub(b1.Mul(c2))
	Dy := a1.Mul(c2).Sub(c1.Mul(a2))

	if !D.Equal(dec.NewFromInt(0)) {
		s := Dx.Div(D)
		t := Dy.Div(D)

		if !(s.LessThanOrEqual(ZERO) && s.LessThanOrEqual(ONE) && t.LessThanOrEqual(ZERO) && t.LessThanOrEqual(ONE)) {
			return Point{ZERO, ZERO}, errors.New("no intersection found")
		} else {
			x := (ONE.Sub(s)).Mul(L1.p1.x).Add(s.Mul(L1.p2.x))
			y := (ONE.Sub(s)).Mul(L1.p1.y).Add(s.Mul(L1.p2.y))
			return (Point{x, y}), nil
		}
	} else {
		return Point{ZERO, ZERO}, errors.New("no intersection found")
	}
}

// func (p *Point) ConnectingEdges()

// LineSegment
type LineSegment struct {
	p1 Point
	p2 Point
}

func (ls *LineSegment) L2() dec.Decimal {
	return ls.p1.x.Sub(ls.p2.x).Pow(TWO).Add((ls.p1.y.Sub(ls.p2.y)).Pow(TWO))
}

func (ls *LineSegment) Mid() Point {
	x := (ls.p1.x.Add(ls.p2.x)).Div(TWO)
	y := (ls.p1.y.Add(ls.p2.y)).Div(TWO)
	return Point{x, y}
}

func (ls *LineSegment) Direction(pt Point) bool {
	direction := (&LineSegment{ls.p1, END}).L2().LessThan((&LineSegment{ls.p2, END}).L2())
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
	n int64
}

func FindNextCoordinate(pt Point, theta dec.Decimal) Point {
	// x = pt.x * cos(t) - pt.y * sin(t)
	// y = pt.x * sin(t) + pt.y * cos(t)
	x := pt.x.Mul(theta.Cos()).Add(pt.y.Mul(theta.Sin()).Neg())
	y := pt.x.Mul(theta.Sin()).Add(pt.y.Mul(theta.Cos()))
	return Point{x, y}
}

func (ps PolygonSolver) CreateRegularPolygon() []Point {
	theta := dec.NewFromInt(2).Div(dec.NewFromInt(ps.n)).Mul(dec.NewFromFloat(math.Pi))
	polygon := make([]Point, 0)
	previous := START

	for i := int64(0); i < ps.n; i++ {
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

	for _, v := range polygonSolver.CreateRegularPolygon() {
		fmt.Println(v.ToString())
	}
}

// Define PolygonSolver
// Define main()

// Create all the polygon nodes
// For each polygon node, create a line segment to every other polygon node, and place all the line segments into a slice
// For each line segment,
