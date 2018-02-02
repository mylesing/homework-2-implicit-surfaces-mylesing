#version 300 es

precision highp float;

in vec4 fs_Pos;

out vec4 out_Col;

uniform vec4 u_Resolution;
uniform mat4 u_ViewMatrix;
uniform mat4 u_ProjMatrix;
uniform vec4 u_CameraPos;
uniform float u_Time;

float sceneSDF(vec3 samplePoint);

//////// CONSTANTS ////////

const int MAX_MARCHING_STEPS = 255;
const float MIN_DIST = 0.0;
const float MAX_DIST = 100.0;
const float EPSILON = 0.01;

/////// PHONG SHADER //////
// estimate normal using SDF gradient
vec3 estimateNormal(vec3 p) {
    return normalize(vec3(
        sceneSDF(vec3(p.x + EPSILON, p.y, p.z)) - sceneSDF(vec3(p.x - EPSILON, p.y, p.z)),
        sceneSDF(vec3(p.x, p.y + EPSILON, p.z)) - sceneSDF(vec3(p.x, p.y - EPSILON, p.z)),
        sceneSDF(vec3(p.x, p.y, p.z + EPSILON)) - sceneSDF(vec3(p.x, p.y, p.z - EPSILON))
    ));
}

vec3 calcNormal( in vec3 pos )
{
    vec2 e = vec2(1.0,-1.0)*0.5773*0.0005;
    return normalize( e.xyy * sceneSDF(pos + e.xyy ) +
					  e.yyx * sceneSDF( pos + e.yyx ) + 
					  e.yxy * sceneSDF( pos + e.yxy ) + 
					  e.xxx * sceneSDF( pos + e.xxx ));
    /*
	vec3 eps = vec3( 0.0005, 0.0, 0.0 );
	vec3 nor = vec3(
	    map(pos+eps.xyy).x - map(pos-eps.xyy).x,
	    map(pos+eps.yxy).x - map(pos-eps.yxy).x,
	    map(pos+eps.yyx).x - map(pos-eps.yyx).x );
	return normalize(nor);
	*/
}

/**
 * Lighting contribution of a single point light source via Phong illumination.
 * 
 * The vec3 returned is the RGB color of the light's contribution.
 *
 * k_a: Ambient color
 * k_d: Diffuse color
 * k_s: Specular color
 * alpha: Shininess coefficient
 * p: position of point being lit
 * eye: the position of the camera
 * lightPos: the position of the light
 * lightIntensity: color/intensity of the light
 *
 * See https://en.wikipedia.org/wiki/Phong_reflection_model#Description
 */
vec3 phongContribForLight(vec3 k_d, vec3 k_s, float alpha, vec3 p, vec3 eye,
                          vec3 lightPos, vec3 lightIntensity) {
    vec3 N = estimateNormal(p);
    vec3 L = normalize(lightPos - p);
    vec3 V = normalize(eye - p);
    vec3 R = normalize(reflect(-L, N));
    
    float dotLN = dot(L, N);
    float dotRV = dot(R, V);
    
    if (dotLN < 0.0) {
        // Light not visible from this point on the surface
        return vec3(0.0, 0.0, 0.0);
    } 
    
    if (dotRV < 0.0) {
        // Light reflection in opposite direction as viewer, apply only diffuse
        // component
        return lightIntensity * (k_d * dotLN);
    }
    return lightIntensity * (k_d * dotLN + k_s * pow(dotRV, alpha));
}

/**
 * Lighting via Phong illumination.
 * 
 * The vec3 returned is the RGB color of that point after lighting is applied.
 * k_a: Ambient color
 * k_d: Diffuse color
 * k_s: Specular color
 * alpha: Shininess coefficient
 * p: position of point being lit
 * eye: the position of the camera
 *
 * See https://en.wikipedia.org/wiki/Phong_reflection_model#Description
 */
vec3 phongIllumination(vec3 k_a, vec3 k_d, vec3 k_s, float alpha, vec3 p, vec3 eye) {
    const vec3 ambientLight = 0.5 * vec3(1.0, 1.0, 1.0);
    vec3 color = ambientLight * k_a;
    
    vec3 lightPos = vec3(4.0 * sin(u_Time * 0.01),
                          2.0,
                          4.0 * cos(u_Time * 0.01));
    vec3 lightIntensity = vec3(0.4, 0.4, 0.4);
    
    // color += phongContribForLight(k_d, k_s, alpha, p, eye,
    //                               lightPos,
    //                               lightIntensity);
      
    return color * 2.0;
}

/////// RAY MARCHER ////////

// SDF sphere
float sphereSDF(vec3 samplePoint, float rad) {
    return length(samplePoint) - rad;
}

// SDF torus
float length2( vec2 p )
{
	return sqrt( p.x*p.x + p.y*p.y );
}

float length8( vec2 p )
{
	p = p*p; p = p*p; p = p*p;
	return pow( p.x + p.y, 1.0/8.0 );
}

float sdTorus( vec3 p, vec2 t )
{
  vec2 q = vec2(length2(p.xz)-t.x,p.y);
  return length8(q)-t.y;
}

// union
float opU(float d1, float d2) {
    return min(d1, d2);
}

// subtract
float opS( float d1, float d2 )
{
    return max(-d1, d2);
}

// intersection
float opI( float d1, float d2 )
{
    return max(d1, d2);
}

// rotate
float rotTor(vec3 p, mat4 m, vec2 pos)
{
    vec4 q = inverse(m) * vec4(p, 1.0);
    return sdTorus(q.xyz, pos);
}

// returns a rotation matrix from inputs
mat4 rotate(float angle, float x, float y, float z) {
    float theta = 3.1415 * (angle / 180.0);
    float c = cos(theta);
    float s = sin(theta);

    if (abs(c) < 0.001) {
        c = 0.0;
    }

    if (abs(s) < 0.001) {
        s = 0.0;
    }

    vec4 col1;
    vec4 col2;
    vec4 col3;
    vec4 col4;

    // Get rotation vector:
    if (x == 1.0) { // Rx
        col1 = vec4(1.0, 0.0, 0.0, 0.0);
        col2 = vec4(0.0, c, s, 0.0);
        col3 = vec4(0.0, -s, c, 0.0);
        col4 = vec4(0.0, 0.0, 0.0, 1.0);
    } else if (y == 1.0) { // Ry
        col1 = vec4(c, 0.0, s, 0.0);
        col2 = vec4(0.0, 1.0, 0.0, 0.0);
        col3 = vec4(-s, 0.0, c, 0.0);
        col4 = vec4(0.0, 0.0, 0.0, 1.0);
    } else if (z == 1.0) { //Rz
        col1 = vec4(c, -s, 0.0, 0.0);
        col2 = vec4(s, c, 0.0, 0.0);
        col3 = vec4(0.0, 0.0, 1.0, 0.0);
        col4 = vec4(0.0, 0.0, 0.0, 1.0);
    }

    return mat4(col1, col2, col3, col4);;
}

// power smooth min (k = 8);
float smin( float a, float b, float k )
{
    float res = exp( -k*a ) + exp( -k*b );
    return -log( res )/k;
}

// blend
float opBlend(float s1, float s2, float k)
{
    return smin(s1, s2, k);
}

float sdPlane( vec3 p )
{
	return p.y;
}

float rotSphere(vec3 p, mat4 m, float f)
{
    vec4 q = inverse(m) * vec4(p, 1.0);
    return sphereSDF(q.xyz, f);
}

// Capsule / Line - signed - exact
float sdCapsule( vec3 p, vec3 a, vec3 b, float r )
{
    vec3 pa = p - a;
    vec3 ba = b - a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0);
    return length( pa - ba*h ) - r;
}

float opCheapBend(vec3 p, vec3 a, vec3 b, float r)
{
    float c = cos(20.0*p.y);
    float s = sin(20.0*p.y);
    mat2  m = mat2(c,-s,s,c);
    vec3  q = vec3(m*p.xy,p.z);
    return sdCapsule(q, a, b, r);
}

float sdEllipsoid( in vec3 p, in vec3 r )
{
    return (length( p/r ) - 1.0) * min(min(r.x,r.y),r.z);
}

float opRep(vec3 p, vec3 c)
{
    vec3 q = mod(p,c)-0.5*c;
    mat4 rot = rotate(p.x * p.y * 50.0 * mod(u_Time * 0.01, 100.0), 0.0, 0.0, 1.0);
    vec4 pos = rot * vec4(q, 1.0);
    return sphereSDF(q + vec3(-cos(u_Time * 0.1), sin(u_Time * 0.1), 0.2 * cos(u_Time * 0.1)), 0.1);
}

float udRoundBox( vec3 p, vec3 b, float r )
{
  return length(max(abs(p)-b,0.0))-r;
}

float opRepFloor(vec3 p, vec3 c)
{
    vec3 q = mod(p,c)-0.5*c;
    return  udRoundBox(q, vec3(0.5, 0.05, 0.5), 0.01);
}

// SDF for the scene
float sceneSDF(vec3 samplePoint) {

    // floor
    float res = opRepFloor(samplePoint - vec3(0.0, -2.0, 0.0), vec3(1.1, 0.0, 1.1));

	// head
	samplePoint = samplePoint - vec3(0.45, 0.25, 0.0);
	mat4 rotH = rotate(45.0 + 10.0 * sin(u_Time * 0.05), 0.0, 1.0, 0.0);
	vec4 pos = rotH * vec4(samplePoint, 1.0);
	samplePoint = pos.xyz;
	vec3 headpos = samplePoint;

	float head = sphereSDF(samplePoint, 0.5); 
	head = opS(sphereSDF(samplePoint, 0.45), head);
	mat4 rotT = rotate(-150.0, 0.0, 0.0, 1.0);
	res = opU(res, opBlend(head, rotTor(samplePoint - vec3(-0.28, -0.42, 0.0), rotT, vec2(0.55, 0.005)), 5.5));

	// body 1
	samplePoint = samplePoint + vec3(0.5, 0.75, 0.0);
	float body1 = sphereSDF(samplePoint, 0.6); // opU(res, sphereSDF(samplePoint, 0.5));
	body1 = opS(sphereSDF(samplePoint, 0.55), body1);
	rotT = rotate(-110.0, 0.0, 0.0, 1.0);
	body1 = opBlend(body1, rotTor(samplePoint - vec3(-0.55, -0.3, 0.0), rotT, vec2(0.55, 0.005)), 5.5);
	res = opU(res, body1);

	// body 2
	samplePoint = samplePoint + vec3(0.95, 0.45, 0.2);
	float body2 = sphereSDF(samplePoint, 0.6); // opU(res, sphereSDF(samplePoint, 0.5));
	body2 = opS(sphereSDF(samplePoint, 0.55), body2);
	rotT = rotate(-80.0, 0.0, 0.0, 1.0);
	body2 = opBlend(body2, rotTor(samplePoint - vec3(-0.55, -0.125, 0.0), rotT, vec2(0.65, 0.005)), 4.5);
	res = opU(res, body2);

	// eye 1
	mat4 rotT2 = rotate(-112.0, 0.0, 0.0, 1.0);
	rotT2 = rotate(-60.0, 0.0, 1.0, 0.0) * rotT2;
	vec3 eyePos = headpos + vec3(-0.35, -0.1, 0.4);
	float eye1 = rotTor(eyePos, rotT2, vec2(0.15, 0.025));
	float eyeBall1 = sphereSDF(eyePos, 0.125);
	res = opU(res, eye1);
	res = opU(res, eyeBall1);

	// eye 2
	mat4 rotT3 = rotate(-30.0, 0.0, 1.0, 0.0) * rotT2;
	vec3 eyePos2 = eyePos + vec3(0.35, -0.1, 0.032);
	float eye2 = rotTor(eyePos2, rotT3, vec2(0.12, 0.025));
	float eyeBall2 = sphereSDF(eyePos2, 0.095);
	res = opU(res, eye2);
	res = opU(res, eyeBall2);

    // eye 3
    mat4 rotT4 = rotate(110.0, 0.0, 1.0, 0.0) * rotT2;
    rotT4 = rotate(15.0, 0.0, 0.0, 1.0) * rotT4;
    vec3 eyePos3 = headpos + vec3(-0.20, -0.08, -0.53);
	float eye3 = rotTor(eyePos3, rotT4, vec2(0.15, 0.025));
	float eyeBall3 = sphereSDF(eyePos3, 0.125);
	res = opU(res, eye3);
	res = opU(res, eyeBall3);

    // mouth
    mat4 rotM = rotate(-30.0, 0.0, 0.0, 1.0);
    vec4 mouthPos = rotM * (vec4(headpos, 1.0) + vec4(-0.35, 0.2, 0.0, 0.0));
    float mouth = sdEllipsoid(mouthPos.xyz, vec3(0.5, 0.08 * (abs(sin(u_Time * 0.05 + 0.5)) + 0.1), 0.6));
    res = opS(mouth, res);

    // tail
    mat4 rotTail = rotate(-50.0 + 5.0 * sin(u_Time * 0.5), 0.0, 0.0, 1.0) 
                    * rotate(30.0 * cos(u_Time * 0.5), 0.0, 1.0, 0.0);
    vec4 body2Pos = rotTail * (vec4(samplePoint, 0.0) - vec4(-1.15, -0.2, 0.0, 0.0));
    float sphere = sphereSDF(samplePoint - vec3(-1.05, -0.15, 0.0), 0.5);
    float tail = sdEllipsoid(body2Pos.xyz + vec3(0.5, 0.65, 0.0), vec3(1.2, 0.3, 0.3));
    res = opU(opBlend(sphere, tail, 3.5), res);

    // antennae
    float ear1 = sdEllipsoid(headpos - vec3(0.0, 0.5, 0.5), vec3(0.05, 0.5, 0.05));
    res = opU(res, ear1);

    float ear2 = sdEllipsoid(headpos - vec3(-0.25, 0.5, -0.15), vec3(0.05, 0.5, 0.05));
    res = opU(res, ear2);

	//#define WORMBOI
	#ifdef WORMBOI
	res = sphereSDF(samplePoint);
	for (float i = 1.0; i < 10.0; ++i) {
		float dir = 0.0;
		float offset = 0.0;
		if (mod(i, 2.0) == 0.0) {
			dir = 1.0;
			offset = i;
		} else {
			dir = -1.0;
			offset = i - 1.0;
		}
		res = opBlend(res, sphereSDF(samplePoint + vec3(dir * offset * abs(sin(u_Time)), 0.0, 0.0)));
	}
	#endif 

    return res;
}



// get shortest distance to surface using ray marching
float shortestDistanceToSurface(vec3 eye, vec3 marchingDirection, float start, float end) {
    float depth = start;
    for (int i = 0; i < MAX_MARCHING_STEPS; i++) {
        float dist = sceneSDF(eye + depth * marchingDirection);
        if (dist < EPSILON) {
			return depth;
        }

        // SPHERE CASTING! :D
        depth += dist;
        if (depth >= end) {
            return end;
        }
    }
    return end;
}           

void main() {
	float x = (fs_Pos.x);
	float y = (fs_Pos.y);

	float aspect = (u_Resolution.x / u_Resolution.y);

	vec4 pos = vec4(x * aspect, y, 0.0, 0.0);

	vec3 eye = vec3(0.0, 0.0, -5.0);
	vec3 dir = normalize(pos.xyz - eye);
    
    float dist = shortestDistanceToSurface(eye, dir, MIN_DIST, MAX_DIST);
    
    if (dist >= MAX_DIST - 2.0 * EPSILON) {
        // Didn't hit anything
        // blue
        vec4 col1 = vec4(0.25, 0.3, 0.8, 1.0);

        // pink
        vec4 col2 = vec4(1.2, 0.77, 0.5, 1.0);

        out_Col = mix(col2, col1, fs_Pos.y);
		return;
    }

	// The closest point on the surface to the eyepoint along the view ray
    vec3 p = eye + dist * dir;

    vec3 color = vec3(1.0, 0.0, 0.0);

	float diffuseTerm = dot(normalize(estimateNormal(p)), normalize(vec3(5.0, 5.0, -5.0)));
        // Avoid negative lighting values
        // diffuseTerm = clamp(diffuseTerm, 0, 1);

    float ambientTerm = 0.2;

    float lightIntensity = diffuseTerm + ambientTerm;   

    // irridescent shader
    vec3 a = vec3(0.8, 0.5, 0.4);
    vec3 b = vec3(0.2, 0.4, 0.2);
    vec3 c = vec3(2.0, 1.0, 1.0);
    vec3 d = vec3(0.00, 0.25, 0.25);

    color = a + b * cos(3.0 * 3.1415 * (c * diffuseTerm + d));
    //out_Col = color;
    
    out_Col = vec4(color * lightIntensity, 1.0);
}
