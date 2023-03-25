package geometry

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestNewPointFromIntersection(t *testing.T) {
	type args struct {
		L1 *LineSegment
		L2 *LineSegment
	}
	tests := []struct {
		name string
		args args
		want *Point
	}{
		{
			name: "Intersect at (0.5,0.5)",
			args: args{
				L1: NewLineSegment(zeroZero, oneOne),
				L2: NewLineSegment(oneZero, zeroOne),
			},
			want: halfHalf,
		},
		{
			name: "Lines are parallel",
			args: args{
				L1: NewLineSegment(zeroOne, oneOne),
				L2: NewLineSegment(zeroZero, oneZero),
			},
			want: nil,
		},
		{
			name: "No intersection found",
			args: args{
				L1: NewLineSegment(zeroOne, oneOne),
				L2: NewLineSegment(zeroZero, halfHalf),
			},
			want: nil,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := NewPointFromIntersection(tt.args.L1, tt.args.L2)
			assert.True(t, tt.want.Equal(got))
		})
	}
}
