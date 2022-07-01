package main

import (
	"testing"
)

func TestGetPointToChannelTotalPoints(t *testing.T) {
	precision := 3
	tolerance := 1e-2

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
