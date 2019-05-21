using System;
using UnityEngine;

namespace SimpleKeplerOrbits
{

	public static class KeplerOrbitUtils
	{

		public const double PI_2 = 6.2831853071796d;
		public const double Deg2Rad = 0.017453292519943d;
		public const double Rad2Deg = 57.295779513082d;
		public const double Epsilon = 1.401298E-45d;

		public static double Acosh(double x)
		{
			if (x < 1.0)
			{
				return 0;
			}
			return Math.Log(x + System.Math.Sqrt(x * x - 1.0));
		}


		public static Vector3d ProjectPointOnPlane(Vector3d point, Vector3d planeNormal)
		{
			return point - planeNormal * DotProduct(point, planeNormal);
		}


		public static Vector3 GetRayPlaneIntersectionPoint(Vector3 pointOnPlane, Vector3 normal, Vector3 rayOrigin, Vector3 rayDirection)
		{
			float dotProd = DotProduct(rayDirection, normal);
			if (Math.Abs(dotProd) < Epsilon)
			{
				return new Vector3();
			}
			Vector3 p = rayOrigin + rayDirection * DotProduct((pointOnPlane - rayOrigin), normal) / dotProd;

			p = p - normal * DotProduct(p - pointOnPlane, normal);
			return p;
		}


		public static Vector3 RotateVectorByAngle(Vector3 v, float angleRad, Vector3 n)
		{
			float cosT = Mathf.Cos(angleRad);
			float sinT = Mathf.Sin(angleRad);
			float oneMinusCos = 1f - cosT;
			//matriz de rotação
			float a11 = oneMinusCos * n.x * n.x + cosT;
			float a12 = oneMinusCos * n.x * n.y - n.z * sinT;
			float a13 = oneMinusCos * n.x * n.z + n.y * sinT;
			float a21 = oneMinusCos * n.x * n.y + n.z * sinT;
			float a22 = oneMinusCos * n.y * n.y + cosT;
			float a23 = oneMinusCos * n.y * n.z - n.x * sinT;
			float a31 = oneMinusCos * n.x * n.z - n.y * sinT;
			float a32 = oneMinusCos * n.y * n.z + n.x * sinT;
			float a33 = oneMinusCos * n.z * n.z + cosT;
			return new Vector3(
				v.x * a11 + v.y * a12 + v.z * a13,
				v.x * a21 + v.y * a22 + v.z * a23,
				v.x * a31 + v.y * a32 + v.z * a33
				);
		}


		public static Vector3d RotateVectorByAngle(Vector3d v, double angleRad, Vector3d n)
		{
			double cosT = Math.Cos(angleRad);
			double sinT = Math.Sin(angleRad);
			double oneMinusCos = 1f - cosT;
			//matriz de rotação
			double a11 = oneMinusCos * n.x * n.x + cosT;
			double a12 = oneMinusCos * n.x * n.y - n.z * sinT;
			double a13 = oneMinusCos * n.x * n.z + n.y * sinT;
			double a21 = oneMinusCos * n.x * n.y + n.z * sinT;
			double a22 = oneMinusCos * n.y * n.y + cosT;
			double a23 = oneMinusCos * n.y * n.z - n.x * sinT;
			double a31 = oneMinusCos * n.x * n.z - n.y * sinT;
			double a32 = oneMinusCos * n.y * n.z + n.x * sinT;
			double a33 = oneMinusCos * n.z * n.z + cosT;
			return new Vector3d(
				v.x * a11 + v.y * a12 + v.z * a13,
				v.x * a21 + v.y * a22 + v.z * a23,
				v.x * a31 + v.y * a32 + v.z * a33
				);
		}

	    //produto de dois vetores
		public static float DotProduct(Vector3 a, Vector3 b)
		{
			return a.x * b.x + a.y * b.y + a.z * b.z;
		}

		public static double DotProduct(Vector3d a, Vector3d b)
		{
			return a.x * b.x + a.y * b.y + a.z * b.z;
		}

		public static Vector3 CrossProduct(Vector3 a, Vector3 b)
		{
			return new Vector3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
		}

		public static Vector3d CrossProduct(Vector3d a, Vector3d b)
		{
			return new Vector3d(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
		}

		//calcular velocidade de órbita circular
		public static Vector3 CalcCircleOrbitVelocity(Vector3 attractorPos, Vector3 bodyPos, double attractorMass, double bodyMass, Vector3 orbitNormal, double gConst)
		{
			Vector3 distanceVector = bodyPos - attractorPos;
			float dist = distanceVector.magnitude;
			double MG = attractorMass * gConst;
			double vScalar = Math.Sqrt(MG / dist);
			return CrossProduct(distanceVector, -orbitNormal).normalized * (float)vScalar;
		}

		public static Vector3d CalcCircleOrbitVelocity(Vector3d attractorPos, Vector3d bodyPos, double attractorMass, double bodyMass, Vector3d orbitNormal, double gConst)
		{
			Vector3d distanceVector = bodyPos - attractorPos;
			double dist = distanceVector.magnitude;
			double MG = attractorMass * gConst;
			double vScalar = Math.Sqrt(MG / dist);
			return CrossProduct(distanceVector, -orbitNormal).normalized * vScalar;
		}

	    //calcular os pontos da reta
		public static Vector3[] CalcOrbitPoints(Vector3 attractorPos, Vector3 bodyPos, double attractorMass, double bodyMass, Vector3 relVelocity, double gConst, int pointsCount)
		{
			if (pointsCount < 3 || pointsCount > 1000000)
			{
				return new Vector3[0];
			}
			Vector3 focusPoint = CalcCenterOfMass(attractorPos, attractorMass, bodyPos, bodyMass);
			Vector3 radiusVector = bodyPos - focusPoint;
			float radiusVectorMagnitude = radiusVector.magnitude;
			Vector3 orbitNormal = CrossProduct(radiusVector, relVelocity);
			double MG = (attractorMass + bodyMass) * gConst;
			Vector3 eccVector = CrossProduct(relVelocity, orbitNormal) / (float)MG - radiusVector / radiusVectorMagnitude;
			double focalParameter = orbitNormal.sqrMagnitude / MG;
			float eccentricity = eccVector.magnitude;
			Vector3 minorAxisNormal = -CrossProduct(orbitNormal, eccVector).normalized;
			Vector3 majorAxisNormal = -CrossProduct(orbitNormal, minorAxisNormal).normalized;
			orbitNormal.Normalize();
			double orbitCompressionRatio;
			double semiMajorAxys;
			double semiMinorAxys;
			Vector3 relFocusPoint;
			Vector3 centerPoint;
			if (eccentricity < 1)
			{
				orbitCompressionRatio = 1 - eccentricity * eccentricity;
				semiMajorAxys = focalParameter / orbitCompressionRatio;
				semiMinorAxys = semiMajorAxys * System.Math.Sqrt(orbitCompressionRatio);
				relFocusPoint = (float)semiMajorAxys * eccVector;
				centerPoint = focusPoint - relFocusPoint;
			}
			else
			{
				orbitCompressionRatio = eccentricity * eccentricity - 1.0;
				semiMajorAxys = focalParameter / orbitCompressionRatio;
				semiMinorAxys = semiMajorAxys * System.Math.Sqrt(orbitCompressionRatio);
				relFocusPoint = -(float)semiMajorAxys * eccVector;
				centerPoint = focusPoint - relFocusPoint;
			}

			Vector3[] points = new Vector3[pointsCount];
			double eccVar = 0f;
			for (int i = 0; i < pointsCount; i++)
			{
				Vector3 result = eccentricity < 1 ?
					new Vector3((float)(System.Math.Sin(eccVar) * semiMinorAxys), -(float)(System.Math.Cos(eccVar) * semiMajorAxys)) :
					new Vector3((float)(System.Math.Sinh(eccVar) * semiMinorAxys), (float)(System.Math.Cosh(eccVar) * semiMajorAxys));
				eccVar += Mathf.PI * 2f / (pointsCount - 1);
				points[i] = minorAxisNormal * result.x + majorAxisNormal * result.y + centerPoint;
			}
			return points;
		}

		//calcular centro de massa
		public static Vector3 CalcCenterOfMass(Vector3 pos1, double mass1, Vector3 pos2, double mass2)
		{
			return ((pos1 * (float)mass1) + (pos2 * (float)mass2)) / (float)(mass1 + mass2);
		}

		public static Vector3d CalcCenterOfMass(Vector3d pos1, double mass1, Vector3d pos2, double mass2)
		{
			return ((pos1 * mass1) + (pos2 * mass2)) / (mass1 + mass2);
		}

		public static double ConvertEccentricToTrueAnomaly(double eccentricAnomaly, double eccentricity)
		{
			if (eccentricity < 1d)
			{
				double cosE = System.Math.Cos(eccentricAnomaly);
				double tAnom = System.Math.Acos((cosE - eccentricity) / (1d - eccentricity * cosE));
				if (eccentricAnomaly > Mathf.PI)
				{
					tAnom = PI_2 - tAnom;
				}
				return tAnom;
			}
			else
			{
				double tAnom = System.Math.Atan2(
					System.Math.Sqrt(eccentricity * eccentricity - 1d) * System.Math.Sinh(eccentricAnomaly),
					eccentricity - System.Math.Cosh(eccentricAnomaly)
				);
				return tAnom;
			}
		}

		public static double ConvertTrueToEccentricAnomaly(double trueAnomaly, double eccentricity)
		{
			if (double.IsNaN(eccentricity) || double.IsInfinity(eccentricity))
			{
				return trueAnomaly;
			}
			trueAnomaly = trueAnomaly % PI_2;
			if (eccentricity < 1d)
			{
				if (trueAnomaly < 0)
				{
					trueAnomaly = trueAnomaly + PI_2;
				}
				double cosT2 = Math.Cos(trueAnomaly);
				double eccAnom = Math.Acos((eccentricity + cosT2) / (1d + eccentricity * cosT2));
				if (trueAnomaly > Math.PI)
				{
					eccAnom = PI_2 - eccAnom;
				}
				return eccAnom;
			}
			else
			{
				double cosT = Math.Cos(trueAnomaly);
				double eccAnom = Acosh((eccentricity + cosT) / (1d + eccentricity * cosT)) * System.Math.Sign(trueAnomaly);
				return eccAnom;
			}
		}

		public static double ConvertMeanToEccentricAnomaly(double meanAnomaly, double eccentricity)
		{
			if (eccentricity < 1)
			{
				return KeplerSolver(meanAnomaly, eccentricity);
			}
			else
			{
				return KeplerSolverHyperbolicCase(meanAnomaly, eccentricity);
			}
		}

		public static double ConvertEccentricToMeanAnomaly(double eccentricAnomaly, double eccentricity)
		{
			if (eccentricity < 1)
			{
				return eccentricAnomaly - eccentricity * System.Math.Sin(eccentricAnomaly);
			}
			else
			{
				return System.Math.Sinh(eccentricAnomaly) * eccentricity - eccentricAnomaly;
			}
		}

		public static double CalcTrueAnomalyForDistance(double distance, double eccentricity, double semiMajorAxis)
		{
			if (eccentricity < 1)
			{
				return Math.Acos((semiMajorAxis * (1d - eccentricity * eccentricity) - distance) / (distance * eccentricity));
			}
			else
			{
				return Math.Acos((semiMajorAxis * (eccentricity * eccentricity - 1d) - distance) / (distance * eccentricity));
			}
		}

		public static double KeplerSolver(double meanAnomaly, double eccentricity)
		{
		
			int iterations = eccentricity < 0.4d ? 3 : 5;
			double e = meanAnomaly;
			double esinE;
			double ecosE;
			double deltaE;
			double n;
			for (int i = 0; i < iterations; i++)
			{
				esinE = eccentricity * System.Math.Sin(e);
				ecosE = eccentricity * System.Math.Cos(e);
				deltaE = e - esinE - meanAnomaly;
				n = 1.0 - ecosE;
				e += -5d * deltaE / (n + System.Math.Sign(n) * System.Math.Sqrt(System.Math.Abs(16d * n * n - 20d * deltaE * esinE)));
			}
			return e;
		}

		public static double KeplerSolverHyperbolicCase(double meanAnomaly, double eccentricity)
		{
			double delta = 1d;
	
			double F = System.Math.Log(2d * System.Math.Abs(meanAnomaly) / eccentricity + 1.8d);
			if (double.IsNaN(F) || double.IsInfinity(F))
			{
				return meanAnomaly;
			}
			while (System.Math.Abs(delta) > 1e-005d)
			{
				delta = (eccentricity * (float)System.Math.Sinh(F) - F - meanAnomaly) / (eccentricity * (float)System.Math.Cosh(F) - 1d);
				if (double.IsNaN(delta) || double.IsInfinity(delta))
				{
					return F;
				}
				F -= delta;
			}
			return F;
		}
	}
}
