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
var polygonVerticesTestCases = []struct {
	name            string
	n               int
	polygonVertices ByXY
	lineSegments    []LineSegment
}{
	{
		name: "Test 60-sided polygon",
		n:    60,
		polygonVertices: ByXY{
			{ONE, ZERO},
			{dec.NewFromFloat(0.994521895368273), dec.NewFromFloat(0.104528463267653)},
			{dec.NewFromFloat(0.978147600733805), dec.NewFromFloat(0.207911690817759)},
			{dec.NewFromFloat(0.951056516295153), dec.NewFromFloat(0.309016994374947)},
			{dec.NewFromFloat(0.913545457642600), dec.NewFromFloat(0.406736643075800)},
			{dec.NewFromFloat(0.866025403784438), dec.NewFromFloat(0.500000000000000)},
			{dec.NewFromFloat(0.809016994374947), dec.NewFromFloat(0.587785252292473)},
			{dec.NewFromFloat(0.743144825477394), dec.NewFromFloat(0.669130606358858)},
			{dec.NewFromFloat(0.669130606358858), dec.NewFromFloat(0.743144825477394)},
			{dec.NewFromFloat(0.587785252292473), dec.NewFromFloat(0.809016994374947)},
			{dec.NewFromFloat(0.500000000000000), dec.NewFromFloat(0.866025403784438)},
			{dec.NewFromFloat(0.406736643075800), dec.NewFromFloat(0.913545457642600)},
			{dec.NewFromFloat(0.309016994374947), dec.NewFromFloat(0.951056516295153)},
			{dec.NewFromFloat(0.207911690817759), dec.NewFromFloat(0.978147600733805)},
			{dec.NewFromFloat(0.104528463267653), dec.NewFromFloat(0.994521895368273)},
			{ZERO, ONE},
			{dec.NewFromFloat(-0.104528463267653), dec.NewFromFloat(0.994521895368273)},
			{dec.NewFromFloat(-0.207911690817759), dec.NewFromFloat(0.978147600733805)},
			{dec.NewFromFloat(-0.309016994374947), dec.NewFromFloat(0.951056516295153)},
			{dec.NewFromFloat(-0.406736643075800), dec.NewFromFloat(0.913545457642600)},
			{dec.NewFromFloat(-0.500000000000000), dec.NewFromFloat(0.866025403784438)},
			{dec.NewFromFloat(-0.587785252292473), dec.NewFromFloat(0.809016994374947)},
			{dec.NewFromFloat(-0.669130606358858), dec.NewFromFloat(0.743144825477394)},
			{dec.NewFromFloat(-0.743144825477394), dec.NewFromFloat(0.669130606358858)},
			{dec.NewFromFloat(-0.809016994374947), dec.NewFromFloat(0.587785252292473)},
			{dec.NewFromFloat(-0.866025403784438), dec.NewFromFloat(0.500000000000000)},
			{dec.NewFromFloat(-0.913545457642600), dec.NewFromFloat(0.406736643075800)},
			{dec.NewFromFloat(-0.951056516295153), dec.NewFromFloat(0.309016994374947)},
			{dec.NewFromFloat(-0.978147600733805), dec.NewFromFloat(0.207911690817759)},
			{dec.NewFromFloat(-0.994521895368273), dec.NewFromFloat(0.104528463267653)},
			{ONE.Neg(), ZERO},
			{dec.NewFromFloat(-0.994521895368273), dec.NewFromFloat(-0.104528463267653)},
			{dec.NewFromFloat(-0.978147600733805), dec.NewFromFloat(-0.207911690817759)},
			{dec.NewFromFloat(-0.951056516295153), dec.NewFromFloat(-0.309016994374947)},
			{dec.NewFromFloat(-0.913545457642600), dec.NewFromFloat(-0.406736643075800)},
			{dec.NewFromFloat(-0.866025403784438), dec.NewFromFloat(-0.500000000000000)},
			{dec.NewFromFloat(-0.809016994374947), dec.NewFromFloat(-0.587785252292473)},
			{dec.NewFromFloat(-0.743144825477394), dec.NewFromFloat(-0.669130606358858)},
			{dec.NewFromFloat(-0.669130606358858), dec.NewFromFloat(-0.743144825477394)},
			{dec.NewFromFloat(-0.587785252292473), dec.NewFromFloat(-0.809016994374947)},
			{dec.NewFromFloat(-0.500000000000000), dec.NewFromFloat(-0.866025403784438)},
			{dec.NewFromFloat(-0.406736643075800), dec.NewFromFloat(-0.913545457642600)},
			{dec.NewFromFloat(-0.309016994374947), dec.NewFromFloat(-0.951056516295153)},
			{dec.NewFromFloat(-0.207911690817759), dec.NewFromFloat(-0.978147600733805)},
			{dec.NewFromFloat(-0.104528463267653), dec.NewFromFloat(-0.994521895368273)},
			{ZERO, ONE.Neg()},
			{dec.NewFromFloat(0.104528463267653), dec.NewFromFloat(-0.994521895368273)},
			{dec.NewFromFloat(0.207911690817759), dec.NewFromFloat(-0.978147600733805)},
			{dec.NewFromFloat(0.309016994374947), dec.NewFromFloat(-0.951056516295153)},
			{dec.NewFromFloat(0.406736643075800), dec.NewFromFloat(-0.913545457642600)},
			{dec.NewFromFloat(0.500000000000000), dec.NewFromFloat(-0.866025403784438)},
			{dec.NewFromFloat(0.587785252292473), dec.NewFromFloat(-0.809016994374947)},
			{dec.NewFromFloat(0.669130606358858), dec.NewFromFloat(-0.743144825477394)},
			{dec.NewFromFloat(0.743144825477394), dec.NewFromFloat(-0.669130606358858)},
			{dec.NewFromFloat(0.809016994374947), dec.NewFromFloat(-0.587785252292473)},
			{dec.NewFromFloat(0.866025403784438), dec.NewFromFloat(-0.500000000000000)},
			{dec.NewFromFloat(0.913545457642600), dec.NewFromFloat(-0.406736643075800)},
			{dec.NewFromFloat(0.951056516295153), dec.NewFromFloat(-0.309016994374947)},
			{dec.NewFromFloat(0.978147600733805), dec.NewFromFloat(-0.207911690817759)},
			{dec.NewFromFloat(0.994521895368273), dec.NewFromFloat(-0.104528463267653)},
		},
		lineSegments: []LineSegment{},
	},
}

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

func NewFromBaseCases(a []BasePolygonCase) []PolygonTestCase {
	var results []PolygonTestCase
	for _, v := range a {
		results = append(results, NewFromBaseCase(v))
	}
	return results
}

func NewFromBaseCase(bpc BasePolygonCase) PolygonTestCase {
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

func GetTestCases() []PolygonTestCase {
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

type StringDecimalMap = map[string]dec.Decimal

func TestCreatePolygonVertices(t *testing.T) {
	// Create a new polygon solver struct for every test case
	for _, test := range polygonTestCases {
		ps := PolygonSolver{test.N}
		got := ps.CreatePolygonVertices(precision, tolerance)
		want := test.PolygonVertices

		// Create map with truncated strings as keys and Decimal struct as values for got
		tGot, tWant := make(map[string]Point), make(map[string]Point)
		for _, p := range got {
			tGot[p.StringFixed(precision)] = p
		}
		for _, p := range want {
			tWant[p.StringFixed(precision)] = p
		}

		// Check maps are equal by:
		// 1. Checking all keys in tWant are found in tGot
		counter := 0
		for _, p := range want {
			s := p.StringFixed(precision)
			if _, ok := tGot[s]; !ok {
				counter += 1
				t.Errorf("s = %v in wanted but not found in got", s)
			}
		}
		// 2. Checking all keys in tWant are found in tGot
		for _, p := range got {
			s := p.StringFixed(precision)
			if _, ok := tWant[s]; !ok {
				counter += 1
				t.Errorf("s = %v in got but not found in wanted", s)
			}
		}
		if counter != 0 {
			t.Errorf("Total of %v errors.", counter)
		}

		// 3. Checking that the lengths of the maps are the same
		isLenEqual := len(got) != len(want)
		if isLenEqual {
			t.Errorf("len(got) = %v, \nlen(want) = %v", len(got), len(want))
		}
	}
}

func TestCreatePolygonVerticesConcurrently(t *testing.T) {
	for _, test := range polygonVerticesTestCases {
		ps := PolygonSolver{test.n}
		got := ps.CreatePolygonVerticesConcurrently(precision, tolerance)
		want := test.polygonVertices

		// Create map with truncated strings as keys and Decimal struct as values for got
		tGot, tWant := make(map[string]Point), make(map[string]Point)
		for _, p := range got {
			tGot[p.StringFixed(precision)] = p
		}
		for _, p := range want {
			tWant[p.StringFixed(precision)] = p
		}

		// Check maps are equal by:
		// 1. Checking all keys in tWant are found in tGot
		counter := 0
		for _, p := range want {
			s := p.StringFixed(precision)
			if _, ok := tGot[s]; !ok {
				counter += 1
				t.Errorf("s = %v in wanted but not found in got", s)
			}
		}
		// 2. Checking all keys in tWant are found in tGot
		for _, p := range got {
			s := p.StringFixed(precision)
			if _, ok := tWant[s]; !ok {
				counter += 1
				t.Errorf("s = %v in got but not found in wanted", s)
			}
		}
		if counter != 0 {
			t.Errorf("Total of %v errors.", counter)
		}

		// 3. Checking that the lengths of the maps are the same
		isLenEqual := len(got) != len(want)
		if isLenEqual {
			t.Errorf("len(got) = %v, \nlen(want) = %v", len(got), len(want))
		}
	}
}

func TestDrawLineSegments(t *testing.T) {

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
