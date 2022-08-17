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
	"strings"

	deque "github.com/gammazero/deque"
	dec "github.com/shopspring/decimal"
	combin "gonum.org/v1/gonum/stat/combin"
)

// Constants
var ZERO dec.Decimal = dec.NewFromInt(0)
var ONE dec.Decimal = dec.NewFromInt(1)
var TWO dec.Decimal = dec.NewFromInt(2)
var END Point = Point{ONE.Neg(), ZERO}
var START Point = Point{ONE, ZERO}

func Sqrt(d dec.Decimal, precision uint) dec.Decimal {
	guess := d.Div(dec.NewFromInt(2))
	if guess.LessThanOrEqual(ONE) {
		guess = ONE
	}
	result := guess
	for i := uint(1); i < precision*3; i *= 2 {
		num := result.Mul(result).Sub(d)
		denom := TWO.Mul(result)
		result = result.Sub(num.Div(denom))
	}
	return result.Truncate(int32(precision))
}

func zeroise(a, tolerance dec.Decimal) dec.Decimal {
	if almostEqual(a, ZERO, tolerance) {
		a = ZERO
	}
	return a
}

func almostEqual(a, b, tolerance dec.Decimal) bool {
	return a.Sub(b).Abs().LessThan(tolerance)
}

func isInBetween(left, x, right dec.Decimal) bool {
	return left.LessThanOrEqual(x) && x.LessThanOrEqual(right)
}

// Point definition
// Defines a Point struct that contains the x and y coordinates of a given point
type Point struct {
	X dec.Decimal `json:"x"`
	Y dec.Decimal `json:"y"`
}

func (p Point) String() string {
	return fmt.Sprintf("(%v,%v)", p.X, p.Y)
}

func (p Point) StringFixed(places uint) string {
	intPlaces := int32(places)
	sX := p.X.StringFixed(intPlaces)
	sY := p.Y.StringFixed(intPlaces)
	return fmt.Sprintf("(%v,%v)", sX, sY)
}

func (p Point) ReducedPrecision(precision uint, tolerance dec.Decimal) Point {
	x := zeroise(p.X.Round(int32(precision)), tolerance)
	y := zeroise(p.Y.Round(int32(precision)), tolerance)
	return Point{x, y}
}

func NewFromPointString(s string) (Point, error) {
	// Use regex to capture the first part and second part
	var pt Point
	reX := regexp.MustCompile(`\((.*),`)
	reY := regexp.MustCompile(`,\s*(.*)\)`)
	x, err := dec.NewFromString(reX.FindStringSubmatch(s)[1])
	if err != nil {
		return pt, err
	}
	y, err := dec.NewFromString(reY.FindStringSubmatch(s)[1])
	if err != nil {
		return pt, err
	}
	return Point{x, y}, nil
}

func NewFromIntersection(L1, L2 LineSegment, precision uint, tolerance dec.Decimal) (Point, error) {
	// Constructor for Point via the intersection of two line segments
	a1 := L1.P2.X.Sub(L1.P1.X)
	b1 := L2.P1.X.Sub(L2.P2.X)
	c1 := L2.P1.X.Sub(L1.P1.X)
	a2 := L1.P2.Y.Sub(L1.P1.Y)
	b2 := L2.P1.Y.Sub(L2.P2.Y)
	c2 := L2.P1.Y.Sub(L1.P1.Y)

	d := a1.Mul(b2).Sub(b1.Mul(a2))
	dX := c1.Mul(b2).Sub(b1.Mul(c2))
	dY := a1.Mul(c2).Sub(c1.Mul(a2))

	var pt Point
	s := dX.Div(d)
	t := dY.Div(d)

	if !almostEqual(d, ZERO, tolerance) && (isInBetween(ZERO, s, ONE) && isInBetween(ZERO, t, ONE)) {
		x := (ONE.Sub(s)).Mul(L1.P1.X).Add(s.Mul(L1.P2.X))
		y := (ONE.Sub(s)).Mul(L1.P1.Y).Add(s.Mul(L1.P2.Y))
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

func (ls LineSegment) L2(precision uint) dec.Decimal {
	a := ls.P1.X.Sub(ls.P2.X).Pow(TWO)
	b := ls.P1.Y.Sub(ls.P2.Y).Pow(TWO)
	return Sqrt(a.Add(b), precision)
}

func (ls LineSegment) Mid() Point {
	x := ls.P1.X.Add(ls.P2.X).Div(TWO)
	y := ls.P1.Y.Add(ls.P2.Y).Div(TWO)
	return Point{x, y}
}

func (ls *LineSegment) Direction(precision uint) bool {
	// Provides the direction of the edge in the graph, by comparing the L2 distances
	// (p1-to-end, p2-to-end).
	p1ToEnd := LineSegment{ls.P1, END}.L2(precision)
	p2ToEnd := LineSegment{ls.P2, END}.L2(precision)
	direction := p1ToEnd.GreaterThan(p2ToEnd)
	return direction
}

// Types
type LineToPoints map[LineSegmentString]map[PointString]Point
type PointString string
type AdjacencyList map[PointString]map[PointString]bool
type Indegrees map[PointString]int
type SimpleAdjacencyList map[int]map[int]bool

// To allow sorting of points
type ByXY []Point

func (a ByXY) Len() int {
	return len(a)
}

func (a ByXY) Swap(i, j int) {
	a[i], a[j] = a[j], a[i]
}

func (a ByXY) Less(i, j int) bool {
	if a[i].X.LessThan(a[j].X) {
		return true
	} else if a[i].X.GreaterThan(a[j].X) {
		return false
	} else {
		return a[i].Y.LessThan(a[j].Y)
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
	N         uint
	Precision uint
	Tolerance dec.Decimal
}

func (ps PolygonSolver) FindNextCoordinate(i uint, pt Point, theta dec.Decimal) Point {
	// Returns the next clockwise point from an existing polygon point
	product := theta.Mul(dec.NewFromInt(int64(i)))
	x := pt.X.Mul(product.Cos()).Sub(pt.Y.Mul(product.Sin()))
	y := pt.X.Mul(product.Sin()).Add(pt.Y.Mul(product.Cos()))
	return Point{x, y}
}

func (ps PolygonSolver) CreatePolygonVertices() ByXY {
	// Returns the vertices of a regular polygon
	theta := TWO.Div(dec.NewFromInt(int64(ps.N))).Mul(dec.NewFromFloat(math.Pi))
	polygon := make(ByXY, 0)

	for i := uint(0); i < ps.N; i++ {
		point := ps.FindNextCoordinate(i, START, theta)
		polygon = append(polygon, point)
	}
	return polygon
}

func (ps PolygonSolver) DrawLineSegments(vertices ByXY) []LineSegment {
	// Returns the list of all basic (unbroken) line segments between all pairs of polygon vertices
	c := combin.Combinations(int(ps.N), 2)
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
		p1str := PointString(line.P1.StringFixed(ps.Precision))
		p2str := PointString(line.P2.StringFixed(ps.Precision))
		result[lineStr][p1str] = line.P1
		result[lineStr][p2str] = line.P2
	}

	// Intersect every pair of line segments
	c := combin.Combinations(len(lineSegments), 2)
	for _, v := range c {
		first, second := lineSegments[v[0]], lineSegments[v[1]]
		point, err := NewFromIntersection(first, second, ps.Precision, ps.Tolerance)

		// Line segments may not intersect
		if err == nil {
			pointStr := PointString(point.StringFixed(ps.Precision))
			result[first.StringFixed(ps.Precision)][pointStr] = point
			result[second.StringFixed(ps.Precision)][pointStr] = point
		}
	}
	return result
}

func (ps PolygonSolver) CreateGraph(lineToPoints LineToPoints) (int, AdjacencyList, Indegrees) {
	// Creates the adjacency list of the graph by:
	// 1. Creating all the smaller line segments that make up the edges of the graph
	// 2. Adding the points to the adjacency list
	fmt.Printf("Creating graph for n=%d...\n", ps.N)

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

			// Standardise first and second
			first := PointString(current.StringFixed(ps.Precision))
			second := PointString(previous.StringFixed(ps.Precision))
			if edge.Direction(ps.Precision) {
				first, second = second, first
			}

			// Add to adjacency list and increment value of indegrees
			adjacencyList[first][second] = true
			indegrees[second] += 1
		}
	}

	// Check that the destination nodes exist as a key in the adjacency list
	// TODO: This check will not be required if we only use CreateSimpleGraph
	nEdges := 0
	for _, pts := range adjacencyList {
		for p := range pts {
			_, ok := adjacencyList[p]
			if !ok {
				err := fmt.Errorf("unmatched node found in adjacency destination but not source")
				panic(err)
			}
		}
		nEdges += len(pts)
	}

	return nEdges, adjacencyList, indegrees
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

func (ps PolygonSolver) GetTheoreticalNodes() int {
	// Returns the theoretical number of nodes in regular n-gon with all diagonals drawn
	// Reference: https://oeis.org/A007569
	n := int(ps.N)

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

func (ps PolygonSolver) CheckNodes(a AdjacencyList) error {
	// Checks that the total number of generated nodes are correct
	total := ps.GetTheoreticalNodes()
	if total != len(a) {
		return fmt.Errorf("expected %d nodes, got %d nodes", total, len(a))
	}
	return nil
}

func (ps PolygonSolver) GetTheoreticalRegions() int {
	// Returns the theoretical number of regions in regular n-gon with all diagonals drawn.
	// Reference: https://oeis.org/A007678
	n := int(ps.N)

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

func (ps PolygonSolver) GetTheoreticalEdges() int {
	// Returns the theoretical number of line segments in regular n-gon with all diagonals drawn
	// Reference: https://oeis.org/A135565
	return ps.GetTheoreticalNodes() + ps.GetTheoreticalRegions() - 1
}

func (ps PolygonSolver) CheckEdges(nEdges int) error {
	// Checks that the total number of generated edges are correct
	total := ps.GetTheoreticalEdges()
	if total != nEdges {
		return fmt.Errorf("expected %d edges, got %d edges", total, nEdges)
	}
	return nil
}

func (ps PolygonSolver) CheckAll(nEdges int, a AdjacencyList) error {
	var err error
	err = ps.CheckNodes(a)
	if err != nil {
		panic(err)
	}
	err = ps.CheckEdges(nEdges)
	if err != nil {
		panic(err)
	}
	return nil
}

func (ps PolygonSolver) CheckNodesSimple(g SimpleAdjacencyList) error {
	// Checks that the total number of generated nodes are correct
	want := ps.GetTheoreticalNodes()
	got := len(g)
	if want != got {
		return fmt.Errorf("expected %d nodes, got %d nodes", want, got)
	}
	return nil
}

func (ps PolygonSolver) CheckEdgesSimple(g SimpleAdjacencyList) error {
	// Checks that the total number of generated edges are correct
	want := ps.GetTheoreticalEdges()
	got := 0
	for _, v := range g {
		got += len(v)
	}

	if want != got {
		return fmt.Errorf("expected %d edges, got %d edges", want, got)
	}
	return nil

}

func (ps PolygonSolver) CheckAllSimple(g SimpleAdjacencyList) error {
	var err error
	err = ps.CheckNodesSimple(g)
	if err != nil {
		panic(err)
	}
	err = ps.CheckEdgesSimple(g)
	if err != nil {
		panic(err)
	}
	return nil
}

func (ps PolygonSolver) CreateSimpleGraph(a AdjacencyList, order []PointString) SimpleAdjacencyList {
	// Returns the simplified adjacency list that labels each node with their topological order
	// from 0 to n-1, 0 being the start point and n-1 being the end point.

	// Map point to their respective index in the topological order list
	pointOrdering := make(map[PointString]int)
	for i, point := range order {
		pointOrdering[point] = i
	}

	graph := make(map[int]map[int]bool)
	for k, v := range a {
		downstream := make(map[int]bool)
		for pointStr := range v {
			downstream[pointOrdering[pointStr]] = true
		}
		graph[pointOrdering[k]] = downstream
	}
	return graph
}

func (ps PolygonSolver) GetTopologicalOrdering(a AdjacencyList, i Indegrees) []PointString {
	// Returns the topological order of the graph
	fmt.Println("Getting topological order...")
	var queue deque.Deque[PointString]

	// Make copies of arguments
	adjacencyListCopy := make(AdjacencyList)
	for k, v := range a {
		adjacencyListCopy[k] = v
	}
	indegreesCopy := make(Indegrees)
	for k, v := range i {
		indegreesCopy[k] = v
	}

	pointStr := PointString(START.StringFixed(ps.Precision))
	queue.PushBack(pointStr)
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

func (ps PolygonSolver) Solve(a AdjacencyList, i Indegrees, order []PointString) map[PointString]big.Int {
	// Solves the polygons problem via bottom-up dynamic programming by calculating the value of nodes
	// in topological order by looking up previously calculated values

	// Initialise the bottom-up dynamic programming map
	dp := make(map[PointString]big.Int)
	pointStr := PointString(START.StringFixed(ps.Precision))
	dp[pointStr] = *big.NewInt(1)

	fmt.Println("Summing over the graph...")
	for _, node := range order {
		for nb := range a[node] {
			a, b := dp[nb], dp[node]
			dp[nb] = *new(big.Int).Add(&a, &b)
		}
	}
	return dp
}

func (ps PolygonSolver) SolveSimple(s SimpleAdjacencyList) []big.Int {
	// Solves the polygons problem via bottom-up dynamic programming by calculating the value of nodes
	// in topological order by looking up previously calculated values

	// Initialise the bottom-up dynamic programming slice
	dp := make([]big.Int, len(s))
	dp[0] = *big.NewInt(1)

	fmt.Println("Summing over the graph...")
	for i := 0; i < len(dp); i++ {
		for v := range s[i] {
			a, b := dp[v], dp[i]
			dp[v] = *new(big.Int).Add(&a, &b)
		}
	}
	return dp
}

func (ps PolygonSolver) Result(dp map[PointString]big.Int) big.Int {
	pointStr := PointString(END.StringFixed(ps.Precision))
	return dp[pointStr]
}

func (ps PolygonSolver) ResultSimple(dp []big.Int) big.Int {
	return dp[len(dp)-1]
}

func (ps PolygonSolver) SaveResult(result big.Int, filename string) {
	// Saves the obtained result to disk.
	// Panics if the current results file has a different value.
	jsonFile, err := os.Open(filename)
	if err != nil {
		fmt.Println(err)
	}
	defer jsonFile.Close()

	byteValue, _ := io.ReadAll(jsonFile)

	parsedItem := make(map[uint]json.Number)
	json.Unmarshal([]byte(byteValue), &parsedItem)

	// Check that the entry is the same
	extracted_result, ok := parsedItem[ps.N]

	if ok {
		if string(extracted_result) != result.String() {
			err = fmt.Errorf("the calculated value for n=%v does not equal the existing value", ps.N)
			panic(err)
		}
		fmt.Println("Calculated result matches existing result.")
	} else {
		parsedItem[ps.N] = json.Number(result.String())
		jsonStr, err := json.MarshalIndent(parsedItem, "", "")
		if err != nil {
			panic(err)
		}
		os.WriteFile(filename, jsonStr, 0644)
	}
}

func (ps PolygonSolver) SaveAdjacencyListKeys(a AdjacencyList, filename string) {
	fmt.Printf("Writing to %v...\n", filename)
	jsonFile, err := os.Open(filename)
	if err != nil {
		fmt.Println(err)
	}
	defer jsonFile.Close()

	// Create slice of keys
	adjacencyListKeys := make([]PointString, 0)
	for k := range a {
		adjacencyListKeys = append(adjacencyListKeys, k)
	}

	// Write keys to disk
	jsonStr, err := json.MarshalIndent(adjacencyListKeys, "", "")
	if err != nil {
		panic(err)
	}
	os.WriteFile(filename, jsonStr, 0644)
}

func Run(p uint, t dec.Decimal, n uint) (big.Int, error) {
	ps := PolygonSolver{n, p, t}
	polygonVertices := ps.CreatePolygonVertices()
	lineSegments := ps.DrawLineSegments(polygonVertices)
	lineToPoints := ps.MapLineSegmentsToPoints(lineSegments)
	nEdges, adjacencyList, indegrees := ps.CreateGraph(lineToPoints)
	order := ps.GetTopologicalOrdering(adjacencyList, indegrees)
	dp := ps.Solve(adjacencyList, indegrees, order)
	result := ps.Result(dp)
	fmt.Println(result.String())
	err := ps.CheckAll(nEdges, adjacencyList)
	if err != nil {
		return *big.NewInt(0), err
	}
	ps.SaveResult(result, "results.json")
	return result, nil
}

func RunSimple(p uint, t dec.Decimal, n uint) (big.Int, error) {
	// Creating the simplified graph runs additional checks to ensure that:
	// 1. The total number of nodes are correct
	// 2. The total number of edges are correct
	// 3. The downstream nodes all exist within the simplified graph
	ps := PolygonSolver{n, p, t}
	polygonVertices := ps.CreatePolygonVertices()
	lineSegments := ps.DrawLineSegments(polygonVertices)
	lineToPoints := ps.MapLineSegmentsToPoints(lineSegments)
	_, adjacencyList, indegrees := ps.CreateGraph(lineToPoints)
	order := ps.GetTopologicalOrdering(adjacencyList, indegrees)
	graph := ps.CreateSimpleGraph(adjacencyList, order)
	dp := ps.SolveSimple(graph)
	result := ps.ResultSimple(dp)
	fmt.Println(result.String())
	err := ps.CheckAllSimple(graph)
	if err != nil {
		return *big.NewInt(0), err
	}
	ps.SaveResult(result, "results.json")
	return result, nil
}

func GridSearch(n uint) {
	// Find a combination of p and t that works for a given value of n
	for p := uint(9); p < 20; p++ {
		t := dec.New(1, -int32(p))
		result, err := Run(p, t, n)
		if err != nil || result.Cmp(big.NewInt(0)) == 0 {
			continue
		}
		// If we reach here, we have found a working combination
		// return p
		fmt.Println(p)
	}
	// return 0
}

func main() {
	p := uint(10)
	t := dec.New(1, -int32(p))

	for i := uint(46); i < 60; i += 2 {
		result, err := Run(p, t, i)
		if err != nil || result.Cmp(big.NewInt(0)) == 0 {
			panic(err)
		}
	}
	// GridSearch(60)
}
