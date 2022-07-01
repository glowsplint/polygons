package main

import (
	"math"
	"testing"

	combin "gonum.org/v1/gonum/stat/combin"
)

func d(x, y int) int {
	if x%y == 0 {
		return 1
	} else {
		return 0
	}
}

func theoreticalTotalPoints(n int) int {
	p := func(x, y int) int {
		return int(math.Pow(float64(x), float64(y)))
	}
	interior := combin.Binomial(n, 4) + (-5*p(n, 3)+45*p(n, 2)-70*n+24)/24*d(n, 2) - (3*n/2)*d(n, 4) + (-45*p(n, 2)+262*n)/6*d(n, 6) + 42*n*d(n, 12) + 60*n*d(n, 18) + 35*n*d(n, 24) - 38*n*d(n, 30) - 82*n*d(n, 42) - 330*n*d(n, 60) - 144*n*d(n, 84) - 96*n*d(n, 90) - 144*n*d(n, 120) - 96*n*d(n, 210)
	return interior + n
}

func TestGetPointToChannelTotalPoints(t *testing.T) {
	for i := 4; i < 60; i += 2 {
		ps := PolygonSolver{i}
		lineToPoints := ps.GetLineToPoints()
		adjacencyList := ps.GetFullAdjacencyList(lineToPoints)
		expected := theoreticalTotalPoints(i)
		actual := len(ps.GetPointToChannel(adjacencyList))
		if actual != expected {
			t.Fatalf(`GetPointToChannel(%v) = %v, expected %v`, i, actual, expected)
		}
	}
}
