package main

import (
	"encoding/json"
	"fmt"
	"io/ioutil"
	"os"
	"testing"

	dec "github.com/shopspring/decimal"
)

var precision int32 = 15
var tolerance dec.Decimal = dec.NewFromFloat(1e-10)

type BasePolygonCase struct {
	Name            string        `json:"name"`
	N               int           `json:"n"`
	PolygonVertices [][]float64   `json:"polygonVertices"`
	LineSegments    [][][]float64 `json:"lineSegments"`
}

type PolygonTestCase struct {
	Name            string
	N               int
	PolygonVertices ByXY
	LineSegments    []LineSegment
}

var polygonTestCases = GetTestCases()

func NewFromBaseCase(bpc BasePolygonCase) PolygonTestCase {
	// Creates a PolygonTestCase from a JSON-marshalled BasePolygonCase struct
	var polygonTestCase PolygonTestCase
	var pv ByXY
	var ls []LineSegment

	// Attach constants
	polygonTestCase.N = bpc.N
	polygonTestCase.Name = bpc.Name

	// Create and attach polygon vertices
	for _, v := range bpc.PolygonVertices {
		pv = append(pv, Point{dec.NewFromFloat(v[0]), dec.NewFromFloat(v[1])})
	}
	polygonTestCase.PolygonVertices = pv

	// Create and attach line segments
	for _, v := range bpc.LineSegments {
		ls = append(ls, LineSegment{
			Point{dec.NewFromFloat(v[0][0]), dec.NewFromFloat(v[0][1])},
			Point{dec.NewFromFloat(v[1][0]), dec.NewFromFloat(v[1][1])},
		})
	}
	polygonTestCase.LineSegments = ls
	return polygonTestCase
}

func NewFromBaseCases(a []BasePolygonCase) []PolygonTestCase {
	// Creates a slice of PolygonTestCases from a slice of JSON-marshalled BasePolygonCase structs
	var results []PolygonTestCase
	for _, v := range a {
		results = append(results, NewFromBaseCase(v))
	}
	return results
}

func GetTestCases() []PolygonTestCase {
	// Returns test cases from test_cases.json
	var basePolygonCases []BasePolygonCase

	// Read test_cases.json
	jsonFile, err := os.Open("test_cases.json")

	if err != nil {
		fmt.Println(err)
	}
	defer jsonFile.Close()

	// Create byte array
	byteValue, _ := ioutil.ReadAll(jsonFile)
	json.Unmarshal(byteValue, &basePolygonCases)

	// Create test cases from base cases
	polygonTestCases := NewFromBaseCases(basePolygonCases)
	return polygonTestCases
}

type StringDecimalMap map[string]dec.Decimal

// TODO: Convert this function to accept a []FixedStringer interface rather than ByXY
// So it can be used by LineSegmentSlice

func AreTwoFixedStringerSlicesEqual(a, b []FixedStringer, t *testing.T) bool {
	// Two []Point are considered to be equal if they contain the same Points
	// Create map with truncated strings as keys and Decimal struct as values for got
	t.Helper()

	tA, tB := make(map[string]FixedStringer), make(map[string]FixedStringer)
	for _, p := range a {
		tA[p.StringFixed(precision)] = p
	}
	for _, p := range b {
		tB[p.StringFixed(precision)] = p
	}

	// Check maps are equal by:
	// 1. Checking all keys in tWant are found in tGot
	counter := 0
	for _, p := range b {
		s := p.StringFixed(precision)
		if _, ok := tA[s]; !ok {
			counter += 1
			t.Errorf("s = %v in wanted but not found in got", s)
		}
	}
	// 2. Checking all keys in tWant are found in tGot
	for _, p := range a {
		s := p.StringFixed(precision)
		if _, ok := tB[s]; !ok {
			counter += 1
			t.Errorf("s = %v in got but not found in wanted", s)
		}
	}
	if counter != 0 {
		t.Errorf("Total of %v errors.", counter)
	}

	// 3. Checking that the lengths of the maps are the same
	isLenEqual := len(a) != len(b)
	if isLenEqual {
		t.Errorf("len(got) = %v, \nlen(want) = %v", len(a), len(b))
	}
	return true
}

func TestCreatePolygonVertices(t *testing.T) {
	// Create a new polygon solver struct for every test case
	for _, testCase := range polygonTestCases {
		ps := PolygonSolver{testCase.N, precision, tolerance}

		got, want := make([]FixedStringer, 0), make([]FixedStringer, 0)
		for _, v := range ps.CreatePolygonVertices() {
			got = append(got, v)
		}
		for _, v := range testCase.PolygonVertices {
			want = append(want, v)
		}
		AreTwoFixedStringerSlicesEqual(got, want, t)
	}
}

func benchmarkCreatePolygonVertices(i int, b *testing.B) {
	for n := 0; n < b.N; n++ {
		ps := PolygonSolver{i, precision, tolerance}
		ps.CreatePolygonVertices()
	}
}

func BenchmarkCreatePolygonVertices4(b *testing.B) {
	benchmarkCreatePolygonVertices(4, b)
}

func BenchmarkCreatePolygonVertices15(b *testing.B) {
	benchmarkCreatePolygonVertices(15, b)
}
func BenchmarkCreatePolygonVertices60(b *testing.B) {
	benchmarkCreatePolygonVertices(60, b)
}

// func TestCreatePolygonVerticesConcurrently(t *testing.T) {
// 	for _, testCase := range polygonTestCases {
// 		ps := PolygonSolver{testCase.N, precision, tolerance}
// 		got := ps.CreatePolygonVerticesConcurrently()
// 		want := testCase.PolygonVertices
// 		AreTwoPointSlicesEqual(got, want, t)
// 	}
// }

// func BenchmarkCreatePolygonVerticesConcurrently(b *testing.B) {
// 	for n := 0; n < b.N; n++ {
// 		ps := PolygonSolver{b.N, precision, tolerance}
// 		ps.CreatePolygonVerticesConcurrently()
// 	}
// }

func TestDrawLineSegments(t *testing.T) {
	for _, testCase := range polygonTestCases {
		ps := PolygonSolver{testCase.N, precision, tolerance}
		pv := ps.CreatePolygonVertices()

		got, want := make([]FixedStringer, 0), make([]FixedStringer, 0)
		for _, v := range ps.DrawLineSegments(pv) {
			got = append(got, v)
		}
		for _, v := range testCase.LineSegments {
			want = append(want, v)
		}
		AreTwoFixedStringerSlicesEqual(got, want, t)
	}
}

// func TestGetPointToChannelTotalPoints(t *testing.T) {
// 	for i := 4; i < 60; i += 2 {
// 		ps := PolygonSolver{i}
// 		lineToPoints := ps.MapLineToPoints(precision, tolerance)
// 		adjacencyList := ps.MapFullAdjacencyList(lineToPoints, precision)
// 		expected := TheoreticalTotalPoints(i)
// 		actual := len(ps.MapPointToChannel(adjacencyList))
// 		if actual != expected {
// 			t.Fatalf(`GetPointToChannel(%v) = %v, expected %v`, i, actual, expected)
// 		}
// 	}
// }
