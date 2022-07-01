package main

import (
	"errors"
	"fmt"
	"math"
	"sort"
	"strconv"

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

type LineToPoints = map[LineSegment]mapset.Set[Point]
type PointToLines = map[Point]mapset.Set[LineSegment]
type AdjacencyList = map[Point]mapset.Set[Point]

func roundFloat(val float64, precision int) float64 {
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

func almostEqual(a, b float64, tolerance float64) bool {
	return math.Abs(a-b) <= tolerance
}

func FromIntersection(L1, L2 LineSegment, precision int, tolerance float64) (Point, error) {
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
			if almostEqual(x, 0, tolerance) {
				x = 0
			}
			if almostEqual(y, 0, tolerance) {
				y = 0
			}
			return (Point{roundFloat(x, precision), roundFloat(y, precision)}), nil
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

func FindNextCoordinate(pt Point, theta float64, precision int, tolerance float64) Point {
	x := pt.x*math.Cos(theta) - pt.y*math.Sin(theta)
	y := pt.x*math.Sin(theta) + pt.y*math.Cos(theta)
	if almostEqual(x, 0, tolerance) {
		x = 0
	}
	if almostEqual(y, 0, tolerance) {
		y = 0
	}
	return Point{roundFloat(x, precision), roundFloat(y, precision)}
}

func (ps PolygonSolver) CreateRegularPolygon(precision int, tolerance float64) []Point {
	// Returns the polygon vertices
	theta := 2.0 / float64(ps.n) * math.Pi
	polygon := make([]Point, 0)
	previous := START

	for i := 0; i < ps.n; i++ {
		polygon = append(polygon, previous)
		previous = FindNextCoordinate(previous, theta, precision, tolerance)
	}
	return polygon
}

func (ps PolygonSolver) LineSegments(precision int, tolerance float64) []LineSegment {
	// For all combinations of polygon vertices, create line segments
	c := combin.Combinations(ps.n, 2)
	p := ps.CreateRegularPolygon(precision, tolerance)
	lineSegments := make([]LineSegment, len(c))
	for i, v := range c {
		start, end := p[v[0]], p[v[1]]
		lineSegments[i] = LineSegment{start, end}
	}
	return lineSegments
}

func (ps PolygonSolver) GetLineToPoints(precision int, tolerance float64) LineToPoints {
	// For all combinations of line segments, create intersection points
	// Returns the map of unbroken line segments to points on the line segments
	// TODO: Make concurrent
	lineToPoints := make(LineToPoints)
	lines := ps.LineSegments(precision, tolerance)

	// Add both known polygon points onto the line
	for _, line := range lines {
		lineToPoints[line] = mapset.NewSet(line.p1, line.p2)
	}

	// Intersect all line segments
	c := combin.Combinations(len(lines), 2)
	for _, v := range c {
		first, second := lines[v[0]], lines[v[1]]
		point, err := FromIntersection(first, second, precision, tolerance)
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

func (ps PolygonSolver) GetPointToChannel(adjacencyList AdjacencyList) map[Point]chan int {
	pointToChannel := make(map[Point]chan int)
	for p := range adjacencyList {
		pointToChannel[p] = make(chan int)
	}
	return pointToChannel
}

func GetPointToChannel(precision int, tolerance float64, n int) map[Point]chan int {
	ps := PolygonSolver{n}
	lineToPoints := ps.GetLineToPoints(precision, tolerance)
	adjacencyList := ps.GetFullAdjacencyList(lineToPoints)
	return ps.GetPointToChannel(adjacencyList)
}

func GetTotalPoints(precision int, tolerance float64, n int) int {
	return len(GetPointToChannel(precision, tolerance, n))
}

func GetSortedPoints(precision int, tolerance float64, n int) ByXY {
	pointToChannel := GetPointToChannel(precision, tolerance, n)
	sortedPoints := make(ByXY, 0)
	for k := range pointToChannel {
		sortedPoints = append(sortedPoints, k)
	}
	sort.Sort(sortedPoints)
	return sortedPoints
}

func d(x, y int) int {
	if x%y == 0 {
		return 1
	} else {
		return 0
	}
}

func TheoreticalTotalPoints(n int) int {
	p := func(x, y int) int {
		return int(math.Pow(float64(x), float64(y)))
	}
	interior := combin.Binomial(n, 4) + (-5*p(n, 3)+45*p(n, 2)-70*n+24)/24*d(n, 2) - (3*n/2)*d(n, 4) + (-45*p(n, 2)+262*n)/6*d(n, 6) + 42*n*d(n, 12) + 60*n*d(n, 18) + 35*n*d(n, 24) - 38*n*d(n, 30) - 82*n*d(n, 42) - 330*n*d(n, 60) - 144*n*d(n, 84) - 96*n*d(n, 90) - 144*n*d(n, 120) - 96*n*d(n, 210)
	return interior + n
}

type Result struct {
	p int
	i int
	t float64
	n int
}

func GridSearch() {
	// Find a precision and tolerance that has highest
	bestCombination := Result{}
	for p := 1; p < 16; p++ {
		for i := 1; i < 16; i++ {
			for n := 4; n < 50; n += 2 {
				t, err := strconv.ParseFloat(fmt.Sprintf(`1e-%v`, i), 64)
				if err != nil {
					continue
				}
				totalPoints := GetTotalPoints(p, t, n)
				fmt.Println(p, i, n, totalPoints, TheoreticalTotalPoints(n))
				if totalPoints != TheoreticalTotalPoints(n) {
					break
				}
				if n < bestCombination.n {
					continue
				}
				bestCombination = Result{p, i, t, n}
			}
		}
	}
	fmt.Println(bestCombination)
}

func main() {
	// GridSearch()
	p := 14
	i := 14
	t, _ := strconv.ParseFloat(fmt.Sprintf(`1e-%v`, i), 64)
	n := 12

	sortedPoints := GetSortedPoints(p, t, n)
	for i, v := range sortedPoints {
		fmt.Println(i, v)
	}
	fmt.Println(TheoreticalTotalPoints(n))
}
