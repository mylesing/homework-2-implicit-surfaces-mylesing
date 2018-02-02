#version 300 es

precision highp float;

in vec4 vs_Pos;

out vec4 fs_Pos;

uniform vec4 u_Resolution;

void main() {
	// TODO: Pass relevant info to fragment
	gl_Position = vs_Pos;
	fs_Pos = vs_Pos;
}
