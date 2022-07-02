package main

import (
	"testing"

	dec "github.com/shopspring/decimal"
)

var precision int32 = 15
var tolerance dec.Decimal = dec.NewFromFloat(1e-10)

func TestCreatePolygonVertices(t *testing.T) {
	// Test cases
	tests := []struct {
		name string
		n    int
		want ByXY
	}{
		{
			name: "Test 4-sided polygon",
			n:    4,
			want: ByXY{
				{ONE, ZERO},
				{ZERO, ONE},
				{ONE.Neg(), ZERO},
				{ZERO, ONE.Neg()},
			},
		},
		{
			name: "Test 60-sided polygon",
			n:    60,
			want: ByXY{
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
		},
	}

	// Run all tests
	for _, test := range tests {
		ps := PolygonSolver{test.n}
		got := ps.CreatePolygonVertices(precision, tolerance)
		want := test.want

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
		for _, p := range want {
			s := p.StringFixed(precision)
			if _, ok := tGot[s]; !ok {
				t.Errorf("s = %v in wanted but not found in got", s)
			}
		}
		// 2. Checking all keys in tWant are found in tGot
		for _, p := range got {
			s := p.StringFixed(precision)
			if _, ok := tWant[s]; !ok {
				t.Errorf("s = %v in got but not found in wanted", s)
			}
		}

		// 3. Checking that the lengths of the maps are the same
		isLenEqual := len(got) != len(want)
		if isLenEqual {
			t.Errorf("len(got) = %v, \nlen(want) = %v", len(got), len(want))
		}
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
