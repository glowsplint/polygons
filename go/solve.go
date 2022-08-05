package main

import (
	"encoding/json"
	"errors"
	"fmt"
	"io"
	"math"
	"math/big"
	"os"
	"regexp"
	"sort"
	"strconv"
	"strings"

	deque "github.com/gammazero/deque"
	combin "gonum.org/v1/gonum/stat/combin"
)

// Constants
var END Point = Point{-1, 0}
var START Point = Point{1, 0}

// Point definition
type Point struct {
	X float64 `json:"x"`
	Y float64 `json:"y"`
}

func roundFloat(val float64, precision uint) float64 {
	ratio := math.Pow(10, float64(precision))
	return math.Round(val*ratio) / ratio
}

func (p Point) String() string {
	return fmt.Sprintf("(%v,%v)", p.X, p.Y)
}

func (p Point) StringFixed(places uint) PointString {
	x := roundFloat(p.X, places)
	y := roundFloat(p.Y, places)
	fX := strconv.FormatFloat(x, 'f', int(places), 64)
	fY := strconv.FormatFloat(y, 'f', int(places), 64)
	return PointString(fmt.Sprintf("(%v,%v)", fX, fY))
}

func (p Point) ReducedPrecision(precision uint, tolerance float64) Point {
	x := zeroise(roundFloat(p.X, precision), tolerance)
	y := zeroise(roundFloat(p.Y, precision), tolerance)
	return Point{x, y}
}

func NewFromPointString(s string) (Point, error) {
	// Use regex to capture the first part and second part
	var pt Point
	reX := regexp.MustCompile(`\((.*),`)
	reY := regexp.MustCompile(`,(.*)\)`)
	x, err := strconv.ParseFloat(reX.FindStringSubmatch(s)[1], 64)
	if err != nil {
		return pt, err
	}
	y, err := strconv.ParseFloat(reY.FindStringSubmatch(s)[1], 64)
	if err != nil {
		return pt, err
	}
	return Point{x, y}, nil
}

func zeroise(a float64, tolerance float64) float64 {
	if almostEqual(a, 0, tolerance) {
		a = 0
	}
	return a
}

func NewFromIntersection(L1, L2 LineSegment, precision uint, tolerance float64) (Point, error) {
	// Constructor for Point via the intersection of two line segments
	a1 := L1.P2.X - L1.P1.X
	b1 := L2.P1.X - L2.P2.X
	c1 := L2.P1.X - L1.P1.X
	a2 := L1.P2.Y - L1.P1.Y
	b2 := L2.P1.Y - L2.P2.Y
	c2 := L2.P1.Y - L1.P1.Y

	d := a1*b2 - b1*a2
	dX := c1*b2 - b1*c2
	dY := a1*c2 - c1*a2

	var pt Point
	s := dX / d
	t := dY / d

	if !almostEqual(d, 0, tolerance) && (isInBetween(0, s, 1) && isInBetween(0, t, 1)) {
		x := (1-s)*L1.P1.X + s*L1.P2.X
		y := (1-s)*L1.P1.Y + s*L1.P2.Y
		x = zeroise(x, tolerance)
		y = zeroise(y, tolerance)
		return Point{x, y}.ReducedPrecision(precision, tolerance), nil
	}

	return pt, errors.New("no intersection found")
}

// LineSegment definition
type LineSegment struct {
	P1 Point `json:"p1"`
	P2 Point `json:"p2"`
}

type LineSegmentString string

func (ls LineSegment) String() string {
	return fmt.Sprintf("(%v,%v)", ls.P1.String(), ls.P2.String())
}

func (ls LineSegment) StringFixed(precision uint) LineSegmentString {
	return LineSegmentString(fmt.Sprintf("(%v,%v)", ls.P1.StringFixed(precision), ls.P2.StringFixed(precision)))
}

func (ls LineSegment) L2() float64 {
	a := math.Pow((ls.P1.X - ls.P2.X), 2)
	b := math.Pow((ls.P1.Y - ls.P2.Y), 2)
	return math.Pow(a+b, 0.5)
}

func (ls LineSegment) Mid() Point {
	x := (ls.P1.X + ls.P2.X) / 2
	y := (ls.P1.Y + ls.P2.Y) / 2
	return Point{x, y}
}

func (ls *LineSegment) Direction() bool {
	// Provides the direction of the edge in the graph, by comparing the L2 distances
	// (p1-to-end, p2-to-end).
	p1ToEnd := LineSegment{ls.P1, END}.L2()
	p2ToEnd := LineSegment{ls.P2, END}.L2()
	direction := p1ToEnd > p2ToEnd
	return direction
}

// Types
type LineToPoints map[LineSegmentString]map[PointString]Point
type PointString string
type LineSegments map[LineSegmentString]bool
type AdjacencyList map[PointString]map[PointString]bool
type Indegrees map[PointString]int

// To allow sorting of points
type ByXY []Point

func (a ByXY) Len() int {
	return len(a)
}

func (a ByXY) Swap(i, j int) {
	a[i], a[j] = a[j], a[i]
}

func (a ByXY) Less(i, j int) bool {
	if a[i].X < a[j].X {
		return true
	} else if a[i].X > a[j].X {
		return false
	} else {
		return a[i].Y < a[j].Y
	}
}

func (a ByXY) String() string {
	results := []string{}
	for _, p := range a {
		result := fmt.Sprintf("%v", p)
		results = append(results, result)
	}
	return strings.Join(results, ",")
}

func (s PointString) String() string {
	return fmt.Sprintf("%v", string(s))
}

func (a AdjacencyList) String() string {
	// Prints only keys
	// Convert PointString to Point, sorts and displays
	sortedPoints := make(ByXY, 0)
	for k := range a {
		p, err := NewFromPointString(k.String())
		if err != nil {
			panic(err)
		}
		sortedPoints = append(sortedPoints, p)
	}
	sort.Stable(sortedPoints)
	results := []string{}
	for _, v := range sortedPoints {
		result := fmt.Sprintf("%v", v)
		results = append(results, result)
	}
	return strings.Join(results, "\n")
}

func almostEqual(a, b, tolerance float64) bool {
	return math.Abs(a-b) < tolerance
}

func isInBetween(left, x, right float64) bool {
	return left <= x && x <= right
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
	N         int
	Precision uint
	Tolerance float64
}

func (ps PolygonSolver) FindNextCoordinate(i int, pt Point, theta float64) Point {
	// Returns the next clockwise point from an existing polygon point
	product := theta * float64(i)
	x := (pt.X * (math.Cos(product))) - (pt.Y * (math.Sin(product)))
	y := (pt.X * (math.Sin(product))) + (pt.Y * (math.Cos(product)))
	return Point{x, y}
}

func (ps PolygonSolver) CreatePolygonVertices() ByXY {
	// Returns the vertices of a regular polygon
	theta := 2 / float64(ps.N) * math.Pi
	polygon := make(ByXY, 0)

	for i := 0; i < ps.N; i++ {
		point := ps.FindNextCoordinate(i, START, theta)
		polygon = append(polygon, point)
	}
	return polygon
}

func (ps PolygonSolver) DrawLineSegments(vertices ByXY) []LineSegment {
	// Returns the list of all basic (unbroken) line segments between all pairs of polygon vertices
	c := combin.Combinations(ps.N, 2)
	lineSegments := make([]LineSegment, len(c))
	for i, v := range c {
		start, end := vertices[v[0]], vertices[v[1]]
		start.X = zeroise(start.X, ps.Tolerance)
		start.Y = zeroise(start.Y, ps.Tolerance)
		lineSegments[i] = LineSegment{start, end}
	}
	return lineSegments
}

func (ps PolygonSolver) MapLineSegmentsToPoints(lineSegments []LineSegment) LineToPoints {
	// Returns the mapping of unbroken line segments to all points on the respective line segment.
	result := make(LineToPoints)

	// Add both known polygon points onto the line
	for _, line := range lineSegments {
		lineStr := line.StringFixed(ps.Precision)
		result[lineStr] = make(map[PointString]Point)
		line.P1 = line.P1.ReducedPrecision(ps.Precision, ps.Tolerance)
		line.P2 = line.P2.ReducedPrecision(ps.Precision, ps.Tolerance)
		result[lineStr][line.P1.StringFixed(ps.Precision)] = line.P1
		result[lineStr][line.P2.StringFixed(ps.Precision)] = line.P2
	}

	// Intersect every pair of line segments
	c := combin.Combinations(len(lineSegments), 2)
	for _, v := range c {
		first, second := lineSegments[v[0]], lineSegments[v[1]]
		point, err := NewFromIntersection(first, second, ps.Precision, ps.Tolerance)

		// Line segments may not intersect
		if err == nil {
			pointStr := point.StringFixed(ps.Precision)
			result[first.StringFixed(ps.Precision)][pointStr] = point
			result[second.StringFixed(ps.Precision)][pointStr] = point
		}
	}
	return result
}

func (ps PolygonSolver) CreateGraph(lineToPoints LineToPoints) (LineSegments, AdjacencyList, Indegrees) {
	// Creates the adjacency list of the graph by:
	// 1. Creating all the smaller line segments that make up the edges of the graph
	// 2. Adding the points to the adjacency list
	fmt.Printf("Creating graph for n=%d...\n", ps.N)

	lineSegments := make(LineSegments)
	adjacencyList := make(AdjacencyList)
	indegrees := make(map[PointString]int)

	for _, points := range lineToPoints {
		// Create base adjacency list with all keys
		for pointStr := range points {
			if _, ok := adjacencyList[pointStr]; !ok {
				adjacencyList[pointStr] = make(map[PointString]bool)
				indegrees[pointStr] = 0
			}
		}

		// Sort points by x and y coordinates
		sortedPoints := make(ByXY, 0)
		for _, point := range points {
			sortedPoints = append(sortedPoints, point)
		}
		sort.Stable(sortedPoints)

		for i := 0; i < len(sortedPoints)-1; i++ {
			previous, current := sortedPoints[i], sortedPoints[i+1]

			// Create line segment for each successive pair of points on the unbroken line segment
			edge := LineSegment{previous, current}
			lineSegments[edge.StringFixed(ps.Precision)] = true

			// Standardise first and second
			first := current.StringFixed(ps.Precision)
			second := previous.StringFixed(ps.Precision)
			if edge.Direction() {
				first, second = second, first
			}

			// Add to adjacency list and increment value of indegrees
			adjacencyList[first][second] = true
			indegrees[second] += 1
		}
	}

	return lineSegments, adjacencyList, indegrees
}

func d(x, y int) int {
	if x%y == 0 {
		return 1
	}
	return 0
}
func p(x, y int) int {
	return int(math.Pow(float64(x), float64(y)))
}

func (ps PolygonSolver) GetTheoreticalTotalPoints() int {
	// Returns the theoretical number of nodes in regular n-gon with all diagonals drawn
	// Reference: https://oeis.org/A007569
	n := ps.N

	// tN denotes the term with divisible check in N
	t2 := (-5*p(n, 3) + 45*p(n, 2) - 70*n + 24) / 24 * d(n, 2)
	t4 := (3 * n / 2) * d(n, 4)
	t6 := (-45*p(n, 2) + 262*n) / 6 * d(n, 6)
	t12 := 42 * n * d(n, 12)
	t18 := 60 * n * d(n, 18)
	t24 := 35 * n * d(n, 24)
	t30 := 38 * n * d(n, 30)
	t42 := 82 * n * d(n, 42)
	t60 := 330 * n * d(n, 60)
	t84 := 144 * n * d(n, 84)
	t90 := 96 * n * d(n, 90)
	t120 := 144 * n * d(n, 120)
	t210 := 96 * n * d(n, 210)
	interior := combin.Binomial(n, 4) + t2 - t4 + t6 + t12 + t18 + t24 - t30 - t42 - t60 - t84 - t90 - t120 - t210
	return interior + n
}

func (ps PolygonSolver) CheckIntersections(adjacencyList AdjacencyList) error {
	// Checks that the total number of generated points are correct
	total := ps.GetTheoreticalTotalPoints()
	if total != len(adjacencyList) {
		return fmt.Errorf("expected %d points, got %d points", total, len(adjacencyList))
	}
	return nil
}

func (ps PolygonSolver) GetTheoreticalRegions() int {
	// Returns the theoretical number of regions in regular n-gon with all diagonals drawn.
	// Reference: https://oeis.org/A007678
	n := ps.N

	// tN denotes the term with divisible check in N
	t0 := (p(n, 4) - 6*p(n, 3) + 23*p(n, 2) - 42*n + 24) / 24
	t2 := (-5*p(n, 3) + 42*p(n, 2) - 40*n - 48) / 48 * d(n, 2)
	t4 := -3 * n / 4 * d(n, 4)
	t6 := (-53*p(n, 2) + 310*n) / 12 * d(n, 6)
	t12 := 49 * n / 2 * d(n, 12)
	t18 := 32 * n * d(n, 18)
	t24 := 19 * n * d(n, 24)
	t30 := -36 * n * d(n, 30)
	t42 := -50 * n * d(n, 42)
	t60 := -190 * n * d(n, 60)
	t84 := -78 * n * d(n, 84)
	t90 := -48 * n * d(n, 90)
	t120 := -78 * n * d(n, 120)
	t210 := -48 * n * d(n, 210)
	return t0 + t2 + t4 + t6 + t12 + t18 + t24 + t30 + t42 + t60 + t84 + t90 + t120 + t210
}

func (ps PolygonSolver) GetTheoreticalLineSegments() int {
	// Returns the theoretical number of line segments in regular n-gon with all diagonals drawn
	// Reference: https://oeis.org/A135565
	return ps.GetTheoreticalTotalPoints() + ps.GetTheoreticalRegions() - 1
}

func (ps PolygonSolver) CheckLineSegments(lineSegments LineSegments) error {
	total := ps.GetTheoreticalLineSegments()
	if total != len(lineSegments) {
		return fmt.Errorf("expected %d line segments, got %d line segments", total, len(lineSegments))
	}
	return nil
}

func (ps PolygonSolver) CheckAll(lineSegments LineSegments, adjacencyList AdjacencyList) error {
	var err error
	err = ps.CheckIntersections(adjacencyList)
	if err != nil {
		panic(err)
	}
	err = ps.CheckLineSegments(lineSegments)
	if err != nil {
		panic(err)
	}
	return nil
}

func (ps PolygonSolver) GetTopologicalOrdering(adjacencyList AdjacencyList, indegrees Indegrees) []PointString {
	// Returns the topological order of the graph
	fmt.Println("Getting topological order...")
	var queue deque.Deque[PointString]

	// Make copies of arguments
	adjacencyListCopy := make(AdjacencyList)
	for k, v := range adjacencyList {
		adjacencyListCopy[k] = v
	}
	indegreesCopy := make(Indegrees)
	for k, v := range indegrees {
		indegreesCopy[k] = v
	}

	queue.PushBack(START.StringFixed(ps.Precision))
	order := make([]PointString, 0)

	for queue.Len() > 0 {
		n := queue.PopFront()

		for nb := range adjacencyListCopy[n] {
			indegreesCopy[nb] -= 1
			if indegreesCopy[nb] == 0 {
				queue.PushBack(nb)
			}
		}

		delete(adjacencyListCopy, n)
		order = append(order, n)
	}

	fmt.Println("Topological ordering complete.")
	return order
}

func (ps PolygonSolver) Solve(adjacencyList AdjacencyList, indegrees Indegrees, order []PointString) map[PointString]big.Int {
	// Solves the polygons problem progressively (via bottom-up dynamic programming) by:
	// 1. Topologically sorting the nodes in the graph
	// 2. Calculating the value of nodes by looking up previously calculated values

	// Initialise the bottom-up dynamic programming map
	dp := make(map[PointString]big.Int)
	for node := range adjacencyList {
		dp[node] = *big.NewInt(0)
	}
	dp[START.StringFixed(ps.Precision)] = *big.NewInt(1)

	fmt.Println("Summing over the graph...")
	for _, node := range order {
		for nb := range adjacencyList[node] {
			a, b := dp[nb], dp[node]
			dp[nb] = *new(big.Int).Add(&a, &b)
		}
	}
	return dp
}

func (ps PolygonSolver) Result(dp map[PointString]big.Int) big.Int {
	return dp[END.StringFixed(ps.Precision)]
}

func (ps PolygonSolver) SaveResult(result big.Int, filename string) {
	// Saves the obtained result to disk.
	// Raises an exception if the current results file has a different value.
	jsonFile, err := os.Open(filename)
	if err != nil {
		fmt.Println(err)
	}
	defer jsonFile.Close()

	byteValue, _ := io.ReadAll(jsonFile)

	parsedItem := make(map[string]string)
	json.Unmarshal([]byte(byteValue), &parsedItem)

	// Check that the entry is the same
	s := strconv.Itoa(ps.N)
	extracted_result, ok := parsedItem[s]

	if !ok {
		parsedItem[s] = result.String()
		jsonStr, err := json.MarshalIndent(parsedItem, "", "")
		if err != nil {
			panic(err)
		}
		os.WriteFile(filename, jsonStr, 0644)
	} else if extracted_result != result.String() {
		err = fmt.Errorf("the calculated value for n=%v does not equal the existing value", ps.N)
		panic(err)
	}
	fmt.Println("Calculated result matches existing result.")
}

func (ps PolygonSolver) SaveAdjacencyListKeys(adjacencyList AdjacencyList, filename string) {
	fmt.Printf("Writing to %v...\n", filename)
	jsonFile, err := os.Open(filename)
	if err != nil {
		fmt.Println(err)
	}
	defer jsonFile.Close()

	// Create slice of keys
	adjacencyListKeys := make([]PointString, 0)
	for k := range adjacencyList {
		adjacencyListKeys = append(adjacencyListKeys, k)
	}

	// Write keys to disk
	jsonStr, err := json.MarshalIndent(adjacencyListKeys, "", "")
	if err != nil {
		panic(err)
	}
	os.WriteFile(filename, jsonStr, 0644)
}

func run(p uint, t float64, n int) (big.Int, error) {
	ps := PolygonSolver{n, p, t}
	polygonVertices := ps.CreatePolygonVertices()
	unbrokenLineSegments := ps.DrawLineSegments(polygonVertices)
	lineToPoints := ps.MapLineSegmentsToPoints(unbrokenLineSegments)
	lineSegments, adjacencyList, indegrees := ps.CreateGraph(lineToPoints)
	order := ps.GetTopologicalOrdering(adjacencyList, indegrees)
	dp := ps.Solve(adjacencyList, indegrees, order)
	result := ps.Result(dp)
	fmt.Println(result.String())
	ps.SaveAdjacencyListKeys(adjacencyList, "adjacencyList.json")
	err := ps.CheckAll(lineSegments, adjacencyList)
	if err != nil {
		return *big.NewInt(0), err
	}
	return result, nil
	// ps.SaveResult(result, "results.json")
}

func GridSearch(n int) (uint, int) {
	// Find a combination of p and t that works for a given value of n
	for p := uint(5); p < 20; p++ {
		for k := 5; k < 20; k++ {
			t := math.Pow10(-k)
			result, err := run(p, t, n)
			if err != nil || result.Cmp(big.NewInt(0)) == 0 {
				continue
			}
			// If we reach here, we have found a working combination
			return p, k
		}
	}
	return 0, 0
}

func main() {
	p, t := uint(6), math.Pow10(-5)
	for i := 54; i < 100; i += 2 {
		result, err := run(p, t, i)
		if err != nil || result.Cmp(big.NewInt(0)) == 0 {
			panic(err)
		}
	}
	// p, k := GridSearch(60)
	// fmt.Println(p, k)
}
