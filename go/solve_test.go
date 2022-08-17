package main

import (

	// "reflect"
	"testing"

	dec "github.com/shopspring/decimal"
)

var precision uint = 15
var tolerance dec.Decimal = dec.New(1, -int32(precision))

type FixedStringer interface {
	StringFixed(precision uint) string
}

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
	// 1. Checking all keys in tA are found in tB
	counter := 0
	for _, p := range b {
		s := p.StringFixed(precision)
		if _, ok := tA[s]; !ok {
			counter += 1
			t.Errorf("s = %v in wanted but not found in got", s)
		}
	}
	// 2. Checking all keys in tB are found in tA
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

func TestNewFromPointString(t *testing.T) {
	testCases := []struct {
		input string
		want  Point
	}{
		{"(1.0, 0.0)", Point{dec.NewFromFloat(1.0), dec.NewFromFloat(0.0)}},
		{"(0.0, 1.0)", Point{dec.NewFromFloat(0.0), dec.NewFromFloat(1.0)}},
		{"(-1.0, 0.0)", Point{dec.NewFromFloat(-1.0), dec.NewFromFloat(0.0)}},
		{"(0.0, -1.0)", Point{dec.NewFromFloat(0.0), dec.NewFromFloat(-1.0)}},
	}
	for _, testCase := range testCases {
		got, err := NewFromPointString(testCase.input)
		if err != nil {
			t.Fatal(err)
		}
		// Points are not directly comparable as the struct definition contains pointers
		// We compare their string equivalents instead
		if got.String() != testCase.want.String() {
			t.Errorf("got = %v, want = %v", got, testCase.want)
		}
	}
}

func TestCreatePolygonVertices(t *testing.T) {
	// Create a new polygon solver struct for every test case
	testCases := []struct {
		N        uint
		Vertices []string
	}{
		{4, []string{
			"(1.0, 0.0)",
			"(0.0, 1.0)",
			"(-1.0, 0.0)",
			"(0.0, -1.0)",
		}},
	}

	for _, testCase := range testCases {
		ps := PolygonSolver{testCase.N, precision, tolerance}

		got, want := make([]FixedStringer, 0), make([]FixedStringer, 0)
		for _, v := range ps.CreatePolygonVertices() {
			got = append(got, v)
		}
		for _, v := range testCase.Vertices {
			point, err := NewFromPointString(v)
			if err != nil {
				t.Fatal(err)
			}
			want = append(want, point)
		}
		AreTwoFixedStringerSlicesEqual(got, want, t)
	}
}

func benchmarkCreatePolygonVertices(i int, b *testing.B) {
	for n := 0; n < b.N; n++ {
		ps := PolygonSolver{uint(i), precision, tolerance}
		ps.CreatePolygonVertices()
	}
}

func BenchmarkCreatePolygonVertices4(b *testing.B) {
	benchmarkCreatePolygonVertices(4, b)
}

func BenchmarkCreatePolygonVertices15(b *testing.B) {
	benchmarkCreatePolygonVertices(15, b)
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

// func TestDrawLineSegments(t *testing.T) {
// 	for _, testCase := range polygonTestCases {
// 		ps := PolygonSolver{testCase.N, precision, tolerance}
// 		pv := testCase.PolygonVertices

// 		got, want := make([]FixedStringer, 0), make([]FixedStringer, 0)
// 		for _, v := range ps.DrawLineSegments(pv) {
// 			got = append(got, v)
// 		}
// 		for _, v := range testCase.LineSegments {
// 			want = append(want, v)
// 		}
// 		AreTwoFixedStringerSlicesEqual(got, want, t)
// 	}
// }

// func TestMapLineToPoints4(t *testing.T) {
// 	for _, testCase := range polygonTestCases {
// 		ps := PolygonSolver{testCase.N, precision, tolerance}
// 		lineSegments := testCase.LineSegments
// 		got := ps.MapLineSegmentsToPoints(lineSegments)
// 		want := testCase.LineToPoints

// 		// TODO: Deep equal will not work because points contain Decimal which contains pointers
// 		// Need to convert them into string and do string comparison
// 		if !reflect.DeepEqual(got, want) && ps.N == 4 {
// 			t.Errorf("got = %v, \nwant = %v", got, want)
// 		}
// 	}
// }

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

// func TestGetTheoreticalTotalPoints(t *testing.T) {
// 	for _, testCase := range polygonTestCases {
// 		ps := PolygonSolver{testCase.N, precision, tolerance}
// 		got, want := ps.GetTheoreticalNodes(), testCase.TotalPoints

// 		if got != want {
// 			t.Errorf("got = %v, \nwant = %v", got, want)
// 		}
// 	}
// }

func TestSqrt(t *testing.T) {
	for _, testCase := range []struct {
		Dec       string
		Result    string
		Precision uint
	}{
		{"2", "1.414213562373095", 15},
		{"3", "1.732050807568877", 15},
		{"5", "2.236067977499789", 15},
		{"0.1", "0.316227766016837", 15},
	} {
		d, _ := dec.NewFromString(testCase.Dec)
		want, _ := dec.NewFromString(testCase.Result)
		got := Sqrt(d, testCase.Precision)
		if got.Cmp(want) != 0 {
			t.Errorf("expect %s, got %s", want.String(), got.String())
		}
	}
}
