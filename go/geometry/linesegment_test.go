package geometry

import (
	"testing"

	dec "github.com/shopspring/decimal"
	"github.com/stretchr/testify/assert"
)

func TestNewLineSegment(t *testing.T) {
	type args struct {
		a *Point
		b *Point
	}
	tests := []struct {
		name string
		args args
		want *LineSegment
	}{
		{
			name: "(0,0) and (1,1)",
			args: args{
				a: zeroZero,
				b: oneOne,
			},
			want: &LineSegment{
				P1: zeroZero,
				P2: oneOne,
			},
		},
		{
			name: "(0,0) and (1,1) - reversed",
			args: args{
				b: zeroZero,
				a: oneOne,
			},
			want: &LineSegment{
				P1: zeroZero,
				P2: oneOne,
			},
		},
		{
			name: "(0,0) and (1,0)",
			args: args{
				a: zeroZero,
				b: oneZero,
			},
			want: &LineSegment{
				P1: zeroZero,
				P2: oneZero,
			},
		},
		{
			name: "(0,0) and (1,0) - reversed",
			args: args{
				b: zeroZero,
				a: oneZero,
			},
			want: &LineSegment{
				P1: zeroZero,
				P2: oneZero,
			},
		},
		{
			name: "(0,0) and (0,1)",
			args: args{
				a: zeroZero,
				b: zeroOne,
			},
			want: &LineSegment{
				P1: zeroZero,
				P2: zeroOne,
			},
		},
		{
			name: "(0,0) and (0,1) - reversed",
			args: args{
				b: zeroZero,
				a: zeroOne,
			},
			want: &LineSegment{
				P1: zeroZero,
				P2: zeroOne,
			},
		},
		{
			name: "(0,0) and (0,0)",
			args: args{
				b: zeroZero,
				a: zeroZero,
			},
			want: nil,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := NewLineSegment(tt.args.a, tt.args.b)
			assert.Equal(t, tt.want, got)
		})
	}
}

func TestLineSegment_SquaredL2(t *testing.T) {
	tests := []struct {
		name    string
		ls      *LineSegment
		want    dec.Decimal
		wantErr assert.ErrorAssertionFunc
	}{
		{
			name:    "(0,0) and (1,1)",
			ls:      NewLineSegment(zeroZero, oneOne),
			want:    two,
			wantErr: assert.NoError,
		},
		{
			name:    "(0,0) and (0,0)",
			ls:      NewLineSegment(zeroZero, zeroZero),
			want:    zero,
			wantErr: assert.Error,
		},
		{
			name:    "(0,0) and (0.5,0.5)",
			ls:      NewLineSegment(zeroZero, halfHalf),
			want:    half,
			wantErr: assert.NoError,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := tt.ls.SquaredL2()
			tt.wantErr(t, err)
			assert.True(t, tt.want.Equal(got))
		})
	}
}

func TestLineSegment_Direction(t *testing.T) {
	tests := []struct {
		name    string
		ls      *LineSegment
		want    bool
		wantErr assert.ErrorAssertionFunc
	}{
		{
			name:    "(1,0) further than (0,0) to End",
			ls:      &LineSegment{zeroZero, oneZero},
			want:    false,
			wantErr: assert.NoError,
		},
		{
			name:    "(1,0) further than (0,0) to End",
			ls:      &LineSegment{oneZero, zeroZero},
			want:    true,
			wantErr: assert.NoError,
		},
		{
			name:    "(0,0) is equidistant",
			ls:      &LineSegment{zeroZero, zeroZero},
			want:    false,
			wantErr: assert.Error,
		},
		{
			name:    "p1 is End",
			ls:      &LineSegment{End, zeroZero},
			want:    false,
			wantErr: assert.Error,
		},
		{
			name:    "p2 is End",
			ls:      &LineSegment{zeroZero, End},
			want:    false,
			wantErr: assert.Error,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := tt.ls.Direction()
			tt.wantErr(t, err)
			assert.Equal(t, tt.want, got)
		})
	}
}
