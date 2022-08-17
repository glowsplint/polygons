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

func benchmarkCreatePolygonVertices(i int, b *testing.B) {
	for n := 0; n < b.N; n++ {
		ps := PolygonSolver{uint(i), precision, tolerance}
		ps.CreateExteriorVertices()
	}
}

func BenchmarkCreatePolygonVertices4(b *testing.B) {
	benchmarkCreatePolygonVertices(4, b)
}

func BenchmarkCreatePolygonVertices15(b *testing.B) {
	benchmarkCreatePolygonVertices(15, b)
}

func areTwoFixedStringerSlicesEqual(a, b []FixedStringer, t *testing.T) bool {
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
		want := testCase.want
		// Points are not directly comparable as the struct definition contains pointers
		// We compare their string equivalents instead
		if got.String() != want.String() {
			t.Errorf("got = %v, want = %v", got, testCase.want)
		}
	}
}

func TestCreateExteriorVertices(t *testing.T) {
	testCases := []struct {
		n        uint
		vertices []string
	}{
		{4, []string{
			"(1.0, 0.0)",
			"(0.0, 1.0)",
			"(-1.0, 0.0)",
			"(0.0, -1.0)",
		}},
	}

	for _, testCase := range testCases {
		ps := PolygonSolver{testCase.n, precision, tolerance}

		got, want := make([]FixedStringer, 0), make([]FixedStringer, 0)
		for _, v := range ps.CreateExteriorVertices() {
			got = append(got, v)
		}
		for _, v := range testCase.vertices {
			point, err := NewFromPointString(v)
			if err != nil {
				t.Fatal(err)
			}
			want = append(want, point)
		}
		areTwoFixedStringerSlicesEqual(got, want, t)
	}
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
	testCases := []struct {
		n            uint
		vertices     []string
		lineSegments [][]string
	}{
		{4,
			[]string{
				"(1.0, 0.0)",
				"(0.0, 1.0)",
				"(-1.0, 0.0)",
				"(0.0, -1.0)",
			},
			[][]string{
				{"(1.0, 0.0)", "(0.0, 1.0)"},
				{"(1.0, 0.0)", "(-1.0, 0.0)"},
				{"(1.0, 0.0)", "(0.0, -1.0)"},
				{"(0.0, 1.0)", "(-1.0, 0.0)"},
				{"(0.0, 1.0)", "(0.0, -1.0)"},
				{"(-1.0, 0.0)", "(0.0, -1.0)"},
			}},
	}

	for _, testCase := range testCases {
		ps := PolygonSolver{testCase.n, precision, tolerance}

		// Create vertices from strings
		var pv ByXY
		for _, v := range testCase.vertices {
			pointStr, err := NewFromPointString(v)
			if err != nil {
				t.Fatal(err)
			}
			pv = append(pv, pointStr)
		}

		// Create line segments from strings
		var ls []LineSegment
		for _, v := range testCase.lineSegments {
			p1, err := NewFromPointString(v[0])
			if err != nil {
				t.Fatal(err)
			}
			p2, err := NewFromPointString(v[1])
			if err != nil {
				t.Fatal(err)
			}
			lineSegment := LineSegment{p1, p2}
			ls = append(ls, lineSegment)
		}

		got, want := make([]FixedStringer, 0), make([]FixedStringer, 0)
		for _, v := range ps.DrawLineSegments(pv) {
			got = append(got, v)
		}
		for _, v := range ls {
			want = append(want, v)
		}
		areTwoFixedStringerSlicesEqual(got, want, t)
	}
}

func TestNewFromLineSegmentString(t *testing.T) {
	testCases := []struct {
		input string
		want  LineSegment
	}{
		{"(1, 0), (0, 1)", LineSegment{
			Point{dec.NewFromFloat(1), dec.NewFromFloat(0)},
			Point{dec.NewFromFloat(0), dec.NewFromFloat(1)},
		}},
		{"(-1, 0), (0, -1)", LineSegment{
			Point{dec.NewFromFloat(-1), dec.NewFromFloat(0)},
			Point{dec.NewFromFloat(0), dec.NewFromFloat(-1)},
		}},
	}
	for _, testCase := range testCases {
		got, err := NewFromLineSegmentString(testCase.input)
		if err != nil {
			t.Fatal(err)
		}
		want := testCase.want
		if got.String() != want.String() {
			t.Errorf("got = %v, want = %v", got, testCase.want)
		}
	}
}

func TestMapLineSegmentsToPoints(t *testing.T) {
	testCases := []struct {
		n            uint
		lineSegments [][]string
		lineToPoints map[string][]string
	}{
		{4,
			[][]string{
				{"(1.0, 0.0)", "(0.0, 1.0)"},
				{"(1.0, 0.0)", "(-1.0, 0.0)"},
				{"(1.0, 0.0)", "(0.0, -1.0)"},
				{"(0.0, 1.0)", "(-1.0, 0.0)"},
				{"(0.0, 1.0)", "(0.0, -1.0)"},
				{"(-1.0, 0.0)", "(0.0, -1.0)"},
			},
			map[string][]string{
				"(1, 0), (0, 1)":   {"(0, 1)", "(1, 0)"},
				"(1, 0), (-1, 0)":  {"(-1, 0)", "(0, 0)", "(1, 0)"},
				"(1, 0), (0, -1)":  {"(0, -1)", "(1, 0)"},
				"(0, 1), (-1, 0)":  {"(-1, 0)", "(0, 1)"},
				"(0, 1), (0, -1)":  {"(0, -1)", "(0, 0)", "(0, 1)"},
				"(-1, 0), (0, -1)": {"(-1, 0)", "(0, -1)"},
			}},
	}
	for _, testCase := range testCases {
		ps := PolygonSolver{testCase.n, precision, tolerance}
		// Create line segments from strings
		var ls []LineSegment
		for _, v := range testCase.lineSegments {
			p1, err := NewFromPointString(v[0])
			if err != nil {
				t.Fatal(err)
			}
			p2, err := NewFromPointString(v[1])
			if err != nil {
				t.Fatal(err)
			}
			lineSegment := LineSegment{p1, p2}
			ls = append(ls, lineSegment)
		}

		// Create line to points mapping from strings
		ltp := make(LineToPoints)
		for k, v := range testCase.lineToPoints {
			lineSegment, err := NewFromLineSegmentString(k)
			if err != nil {
				t.Fatal(err)
			}
			key := LineSegmentString(lineSegment.StringFixed(ps.Precision))
			value := make(map[PointString]Point)
			for _, s := range v {
				value[PointString(s)], err = NewFromPointString(s)
				if err != nil {
					t.Fatal(err)
				}
			}
			ltp[key] = value
		}

		got, want := make([]FixedStringer, 0), make([]FixedStringer, 0)
		for _, v := range ps.MapLineSegmentsToPoints(ls) {
			for _, point := range v {
				got = append(got, point)
			}
		}
		for _, v := range ltp {
			for _, point := range v {
				want = append(want, point)
			}
		}
		areTwoFixedStringerSlicesEqual(got, want, t)
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
