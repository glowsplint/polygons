package main

import (
	"errors"
	"fmt"
	"math"
	"os"
	"strings"

	dec "github.com/shopspring/decimal"
	combin "gonum.org/v1/gonum/stat/combin"
)

// Constants
var ZERO dec.Decimal = dec.NewFromInt(0)
var ONE dec.Decimal = dec.NewFromInt(1)
var TWO dec.Decimal = dec.NewFromInt(2)
var END Point = Point{ONE.Neg(), ZERO}
var START Point = Point{ONE, ZERO}

type Point struct {
	x dec.Decimal
	y dec.Decimal
}

func (p Point) String() string {
	return fmt.Sprintf("(%v,%v)", p.x, p.y)
}

func (p Point) StringFixed(places int32) string {
	return fmt.Sprintf("(%v,%v)", p.x.StringFixed(places), p.y.StringFixed(places))
}

type LineToPoints = map[LineSegment]map[string]Point

// type PointToLines = map[Point]mapset.Set[LineSegment]
// type AdjacencyList = map[Point]mapset.Set[Point]

// To allow sorting of points
type ByXY []Point

func (a ByXY) Len() int {
	return len(a)
}

func (a ByXY) Swap(i, j int) {
	a[i], a[j] = a[j], a[i]
}

func (a ByXY) Less(i, j int) bool {
	if a[i].x.LessThan(a[j].x) {
		return true
	} else if a[i].x.GreaterThan(a[j].x) {
		return false
	} else {
		return a[i].y.LessThan(a[j].y)
	}
}

func (a ByXY) String() string {
	results := []string{}
	for _, p := range a {
		result := fmt.Sprintf("%v", p)
		results = append(results, result)
	}
	return strings.Join(results, "\n")
}

func almostEqual(a, b, tolerance dec.Decimal) bool {
	return a.Sub(b).Abs().LessThan(tolerance)
}

func isInBetween(left, x, right dec.Decimal) bool {
	return left.LessThanOrEqual(x) && x.LessThanOrEqual(right)
}

type LineSegment struct {
	p1 Point
	p2 Point
}

func (ls LineSegment) String() string {
	return fmt.Sprintf("(%v,%v)", ls.p1.String(), ls.p2.String())
}

func (ls LineSegment) L2() dec.Decimal {
	a := ls.p1.x.Sub(ls.p2.x).Pow(TWO)
	b := ls.p1.y.Sub(ls.p2.y).Pow(TWO)
	return a.Add(b).Pow(ONE.Div(TWO))
}

func (ls LineSegment) Mid() Point {
	x := ls.p1.x.Add(ls.p2.x).Div(TWO)
	y := ls.p1.y.Add(ls.p2.y).Div(TWO)
	return Point{x, y}
}

func (ls *LineSegment) Direction(pt Point) bool {
	// pt must either be p1 or p2, otherwise the return value does not make sense.
	p1ToEnd := (&LineSegment{ls.p1, END}).L2()
	p2ToEnd := (&LineSegment{ls.p2, END}).L2()
	direction := p1ToEnd.GreaterThan(p2ToEnd)
	if pt == ls.p2 {
		return direction
	} else if pt == ls.p1 {
		return !direction
	}
	panic("The value of argument pt must either be p1 or p2.")
}

func (ls *LineSegment) OtherEnd(pt Point) Point {
	if pt == ls.p1 {
		return ls.p2
	} else if pt == ls.p2 {
		return ls.p1
	}
	panic("The value of argument pt must either be p1 or p2.")
}

func FromIntersection(L1, L2 LineSegment, precision int32, tolerance dec.Decimal) (Point, error) {
	// Constructor for Point via the intersection of two line segments
	a1 := L1.p2.x.Sub(L1.p1.x)
	b1 := L2.p1.x.Sub(L2.p2.x)
	c1 := L2.p1.x.Sub(L1.p1.x)
	a2 := L1.p2.y.Sub(L1.p1.y)
	b2 := L2.p1.y.Sub(L2.p2.y)
	c2 := L2.p1.y.Sub(L1.p1.y)

	D := a1.Mul(b2).Sub(b1.Mul(a2))
	Dx := c1.Mul(b2).Sub(b1.Mul(c2))
	Dy := a1.Mul(c2).Sub(c1.Mul(a2))

	var pt Point
	if !almostEqual(D, ZERO, tolerance) {
		s := Dx.Div(D)
		t := Dy.Div(D)

		if !(isInBetween(ZERO, s, ONE) && isInBetween(ZERO, t, ONE)) {
			return pt, errors.New("no intersection found")
		} else {
			x := ONE.Sub(s).Mul(L1.p1.x).Add(s.Mul(L1.p2.x))
			y := ONE.Sub(s).Mul(L1.p1.y).Add(s.Mul(L1.p2.y))
			if almostEqual(x, ZERO, tolerance) {
				x = ZERO
			}
			if almostEqual(y, ZERO, tolerance) {
				y = ZERO
			}
			return (Point{x.Round(precision), y.Round(precision)}), nil
		}
	} else {
		return pt, errors.New("no intersection found")
	}
}

type LineSegmentSlice []LineSegment

func (lsSlice LineSegmentSlice) String() string {
	results := []string{}
	for _, ls := range lsSlice {
		result := fmt.Sprintf("%v", ls)
		results = append(results, result)
	}
	return strings.Join(results, "\n")
}

type PolygonSolver struct {
	n         int
	precision int32
	tolerance dec.Decimal
}

func (ps PolygonSolver) FindNextCoordinate(pt Point, theta dec.Decimal) Point {
	// Creates the next polygon point from an existing polygon point
	x := pt.x.Mul(theta.Cos()).Sub(pt.y.Mul(theta.Sin()))
	y := pt.x.Mul(theta.Sin()).Add(pt.y.Mul(theta.Cos()))
	if almostEqual(x, ZERO, ps.tolerance) {
		x = ZERO
	}
	if almostEqual(y, ZERO, ps.tolerance) {
		y = ZERO
	}
	return Point{x.Round(ps.precision), y.Round(ps.precision)}
}

func (ps PolygonSolver) CreatePolygonVertices() ByXY {
	// Returns the polygon vertices
	// TODO: Compare speed with a concurrent version by multiplying theta by a constant in a loop
	theta := TWO.Div(dec.NewFromInt(int64(ps.n))).Mul(dec.NewFromFloat(math.Pi))
	polygon := make(ByXY, 0)
	previous := START

	for i := 0; i < ps.n; i++ {
		polygon = append(polygon, previous)
		previous = ps.FindNextCoordinate(previous, theta)
	}
	return polygon
}

func (ps PolygonSolver) GetPolygonVertex(i int, theta dec.Decimal, points chan Point) {
	// Creates a polygon point from the base polygon point
	decI := dec.NewFromInt(int64(i))
	product := decI.Mul(theta)
	x := START.x.Mul(product.Cos()).Sub(START.y.Mul(product.Sin()))
	y := START.x.Mul(product.Sin()).Add(START.y.Mul(product.Cos()))
	if almostEqual(x, ZERO, ps.tolerance) {
		x = ZERO
	}
	if almostEqual(y, ZERO, ps.tolerance) {
		y = ZERO
	}
	points <- Point{x.Round(ps.precision), y.Round(ps.precision)}
}

func (ps PolygonSolver) CreatePolygonVerticesConcurrently() ByXY {
	// Returns the polygon vertices generated concurrently
	theta := TWO.Div(dec.NewFromInt(int64(ps.n))).Mul(dec.NewFromFloat(math.Pi))
	points := make(chan Point)
	polygon := make([]Point, 0)
	for i := 0; i < ps.n; i++ {
		go ps.GetPolygonVertex(i, theta, points)
	}
	for i := 0; i < ps.n; i++ {
		point := <-points
		polygon = append(polygon, point)
	}
	return polygon
}

func (ps PolygonSolver) DrawLineSegments(polygonVertices ByXY) LineSegmentSlice {
	// For all combinations of polygon vertices, create line segments
	c := combin.Combinations(ps.n, 2)
	lineSegments := make([]LineSegment, len(c))
	for i, v := range c {
		start, end := polygonVertices[v[0]], polygonVertices[v[1]]
		lineSegments[i] = LineSegment{start, end}
	}
	return lineSegments
}

// func (ps PolygonSolver) MapLineToPoints(precision int32, tolerance dec.Decimal) LineToPoints {
// 	// For all combinations of line segments, create intersection points
// 	// Returns the map of unbroken line segments to points on the line segments
// 	// TODO: Make concurrent since this is an embarrassingly parallel problem
// 	// TODO: Need to rewrite some part of this function and earlier
// 	lineToPoints := make(LineToPoints)
// 	lines := ps.DrawLineSegments(precision, tolerance)

// 	// Add both known polygon points onto the line
// 	for _, line := range lines {
// 		lineToPoints[line] = mapset.NewSet(line.p1, line.p2)
// 	}

// 	// Intersect all line segments
// 	c := combin.Combinations(len(lines), 2)
// 	for _, v := range c {
// 		first, second := lines[v[0]], lines[v[1]]
// 		point, err := FromIntersection(first, second, precision, tolerance)
// 		if err == nil {
// 			lineToPoints[first].Add(point)
// 			lineToPoints[second].Add(point)
// 		}
// 	}
// 	// TODO: There are more points than there should be, some of them repeated
// 	// Try with .String()
// 	fmt.Println(lineToPoints)
// 	fmt.Println(len(lineToPoints))
// 	return lineToPoints
// }

// func (ps PolygonSolver) MapFullAdjacencyList(lineToPoints LineToPoints, precision int32) AdjacencyList {
// 	// Returns the full adjacency list of nodes and their connected nodes (ignoring direction)
// 	adjacencyList := make(AdjacencyList)
// 	for _, points := range lineToPoints {
// 		// Sort all the points on the unbroken line and get each broken line segment
// 		sortedPoints := points.ToSlice()
// 		sort.Stable(ByXY(sortedPoints))

// 		for i := 1; i < len(sortedPoints); i++ {
// 			// This part is can cause issues because the values are different in decimal form,
// 			// even though the differences are tiny.
// 			// A better way is to use a tree data structure to allow O(log n) searching
// 			// If we truncate it to certain number of digits, we can use O(1) searching
// 			currPt := Point{sortedPoints[i].x.Truncate(precision), sortedPoints[i].y.Truncate(precision)}
// 			prevPt := Point{sortedPoints[i-1].x.Truncate(precision), sortedPoints[i-1].y.Truncate(precision)}
// 			// fmt.Println(currPt, currPt.x.Equal(ZERO), currPt.y.Equal(ONE.Neg()))
// 			// Can confirm that similar value decimals are being placed into the hash
// 			// Next step: check that struct values that comprise decimal
// 			// Unable to do this because struct fields are not exported

// 			if _, ok := adjacencyList[currPt]; ok {
// 				adjacencyList[currPt].Add(prevPt)
// 			} else {
// 				adjacencyList[currPt] = mapset.NewSet(prevPt)
// 			}

// 			if _, ok := adjacencyList[prevPt]; ok {
// 				adjacencyList[prevPt].Add(currPt)
// 			} else {
// 				adjacencyList[prevPt] = mapset.NewSet(currPt)
// 			}
// 		}
// 	}
// 	return adjacencyList
// }

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

// func (ps PolygonSolver) MapPointToChannel(adjacencyList AdjacencyList) map[Point]chan int {
// 	pointToChannel := make(map[Point]chan int)
// 	for p := range adjacencyList {
// 		pointToChannel[p] = make(chan int)
// 	}
// 	return pointToChannel
// }

// func GetLineToPoints(precision int32, tolerance dec.Decimal, n int) LineToPoints {
// 	ps := PolygonSolver{n}
// 	return ps.MapLineToPoints(precision, tolerance)
// }

// func GetFullAdjacencyList(precision int32, tolerance dec.Decimal, n int) AdjacencyList {
// 	ps := PolygonSolver{n}
// 	lineToPoints := ps.MapLineToPoints(precision, tolerance)
// 	return ps.MapFullAdjacencyList(lineToPoints, precision)
// }

// func GetPointToChannel(precision int32, tolerance dec.Decimal, n int) map[Point]chan int {
// 	ps := PolygonSolver{n}
// 	lineToPoints := ps.MapLineToPoints(precision, tolerance)
// 	adjacencyList := ps.MapFullAdjacencyList(lineToPoints, precision)
// 	return ps.MapPointToChannel(adjacencyList)
// }

// func GetTotalPoints(precision int32, tolerance dec.Decimal, n int) int {
// 	return len(GetPointToChannel(precision, tolerance, n))
// }

// func GetSortedPoints(precision int32, tolerance dec.Decimal, n int) ByXY {
// 	pointToChannel := GetPointToChannel(precision, tolerance, n)
// 	sortedPoints := make(ByXY, 0)
// 	for k := range pointToChannel {
// 		sortedPoints = append(sortedPoints, k)
// 	}
// 	sort.Sort(sortedPoints)
// 	return sortedPoints
// }

// func d(x, y int) int {
// 	if x%y == 0 {
// 		return 1
// 	}
// 	return 0
// }

// func TheoreticalTotalPoints(n int) int {
// 	// Returns the theoretical number of total points for the intersections of all line segments
// 	// created by the vertices of a polygon given by the sum of polygon vertices and interior
// 	// intersection points
// 	p := func(x, y int) int {
// 		return int(math.Pow(float64(x), float64(y)))
// 	}

// 	// tN denotes the term with divisible check in N
// 	t2 := (-5*p(n, 3) + 45*p(n, 2) - 70*n + 24) / 24 * d(n, 2)
// 	t4 := (3 * n / 2) * d(n, 4)
// 	t6 := (-45*p(n, 2) + 262*n) / 6 * d(n, 6)
// 	t12 := 42 * n * d(n, 12)
// 	t18 := 60 * n * d(n, 18)
// 	t24 := 35 * n * d(n, 24)
// 	t30 := 38 * n * d(n, 30)
// 	t42 := 82 * n * d(n, 42)
// 	t60 := 330 * n * d(n, 60)
// 	t84 := 144 * n * d(n, 84)
// 	t90 := 96 * n * d(n, 90)
// 	t120 := 144 * n * d(n, 120)
// 	t210 := 96 * n * d(n, 210)
// 	interior := combin.Binomial(n, 4) + t2 - t4 + t6 + t12 + t18 + t24 - t30 - t42 - t60 - t84 - t90 - t120 - t210
// 	return interior + n
// }

// type Result struct {
// 	p int
// 	i int
// 	t dec.Decimal
// 	n int
// }

// func GridSearch() {
// 	// Find a precision and tolerance that has highest
// 	bestCombination := Result{}
// 	for p := 1; p < 16; p++ {
// 		for i := 1; i < 16; i++ {
// 			for n := 4; n < 50; n += 2 {
// 				t := dec.NewFromFloatWithExponent(0.1, int32(-i))
// 				totalPoints := GetTotalPoints(int32(p), t, n)
// 				fmt.Println(p, i, n, totalPoints, TheoreticalTotalPoints(n))
// 				if totalPoints != TheoreticalTotalPoints(n) {
// 					break
// 				}
// 				if n < bestCombination.n {
// 					continue
// 				}
// 				bestCombination = Result{p, i, t, n}
// 			}
// 		}
// 	}
// 	fmt.Println(bestCombination)
// }

func main() {
	// GridSearch()
	p, t := int32(15), dec.NewFromFloat(1e-10)

	ps := PolygonSolver{60, p, t}
	polygonVertices := ps.CreatePolygonVertices()
	lineSegments := ps.DrawLineSegments(polygonVertices)
	file, fileErr := os.Create("results")
	if fileErr != nil {
		fmt.Println(fileErr)
		return
	}
	fmt.Fprintf(file, "%v\n", lineSegments)
}
