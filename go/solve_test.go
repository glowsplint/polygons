package main

import (
	"testing"

	dec "github.com/shopspring/decimal"
)

func TestGetPointToChannelTotalPoints(t *testing.T) {
	precision := int32(3)
	tolerance := dec.NewFromFloatWithExponent(0.1, -2)

	for i := 4; i < 60; i += 2 {
		ps := PolygonSolver{i}
		lineToPoints := ps.GetLineToPoints(precision, tolerance)
		adjacencyList := ps.GetFullAdjacencyList(lineToPoints)
		expected := TheoreticalTotalPoints(i)
		actual := len(ps.GetPointToChannel(adjacencyList))
		if actual != expected {
			t.Fatalf(`GetPointToChannel(%v) = %v, expected %v`, i, actual, expected)
		}
	}
}
