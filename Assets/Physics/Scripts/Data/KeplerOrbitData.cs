﻿#region Copyright
/// Copyright © 2017-2018 Vlad Kirpichenko
/// 
/// Author: Vlad Kirpichenko 'itanksp@gmail.com'
/// Licensed under the MIT License.
/// License: http://opensource.org/licenses/MIT
#endregion

using System;
using UnityEngine;

namespace SimpleKeplerOrbits
{
	/// <summary>
	/// Orbit data container.
	/// Also contains methods for altering and updating orbit state.
	/// </summary>
	[Serializable]
	public class KeplerOrbitData
	{
		public double GravConst = 1;

		/// <summary>
		/// Noraml of ecliptic plane.
		/// </summary>
		public Vector3d EclipticNormal = new Vector3d(0, 0, 1);

		/// <summary>
		/// Up direction on ecliptic plane (y-axis on xy ecliptic plane).
		/// </summary>
		public Vector3d EclipticUp = new Vector3d(0, 1, 0);

		/// <summary>
		/// Right vector on ecliptic plane (x-axis on xy ecliptic plane).
		/// </summary>
		public Vector3d EclipticRight = new Vector3d(1, 0, 0);

		/// <summary>
		/// Body position relatve to attractor or Focal Position.
		/// </summary>
		/// <remarks>
		/// Attractor (focus) is local center of orbit system.
		/// </remarks>
		public Vector3d Position;

		/// <summary>
		/// Magnitude of body position vector.
		/// </summary>
		public double AttractorDistance;

		/// <summary>
		/// Attractor point mass.
		/// </summary>
		public double AttractorMass;

		/// <summary>
		/// Body velocity vector relative to attractor.
		/// </summary>
		public Vector3d Velocity;

		/// <summary>
		/// Gravitational parameter of system.
		/// </summary>
		public double MG;

		public double SemiMinorAxis;
		public double SemiMajorAxis;
		public double FocalParameter;
		public double Eccentricity;
		public double EnergyTotal;
		public double Period;
		public double TrueAnomaly;
		public double MeanAnomaly;
		public double EccentricAnomaly;
		public double SquaresConstant;
		public Vector3d Periapsis;
		public double PeriapsisDistance;
		public Vector3d Apoapsis;
		public double ApoapsisDistance;
		public Vector3d CenterPoint;
		public double OrbitCompressionRatio;
		public Vector3d OrbitNormal;
		public Vector3d SemiMinorAxisBasis;
		public Vector3d SemiMajorAxisBasis;

		/// <summary>
		/// if > 0, then orbit motion is clockwise
		/// </summary>
		public double OrbitNormalDotEclipticNormal;

		/// <summary>
		/// The orbit inclination in radians relative to ecliptic plane.
		/// </summary>
		public double Inclination
		{
			get
			{
				var dot = KeplerOrbitUtils.DotProduct(OrbitNormal, EclipticNormal);
				return Math.Acos(dot);
			}
		}

		/// <summary>
		/// Ascending node longitude in radians.
		/// </summary>
		public double AscendingNodeLongitude
		{
			get
			{
				var ascNodeDir = KeplerOrbitUtils.CrossProduct(EclipticNormal, OrbitNormal).normalized;
				
				var dot = KeplerOrbitUtils.DotProduct(ascNodeDir, EclipticRight);
				if (dot < KeplerOrbitUtils.Epsilon && dot > -KeplerOrbitUtils.Epsilon)
				{
					return 0;
				}
				var angle = Math.Acos(dot);
				if (KeplerOrbitUtils.DotProduct(KeplerOrbitUtils.CrossProduct(ascNodeDir, EclipticRight), EclipticNormal) >= 0)
				{
					angle = KeplerOrbitUtils.PI_2 - angle;
				}
				return angle;
			}
		}

		/// <summary>
		/// Angle between main orbit axis and ecliptic 0 axis in radians.
		/// </summary>
		public double ArgumentOfPerifocus
		{
			get
			{
				var ascNodeDir = KeplerOrbitUtils.CrossProduct(EclipticNormal, OrbitNormal).normalized;
				var dot = KeplerOrbitUtils.DotProduct(ascNodeDir, SemiMajorAxisBasis.normalized);
				if (dot < KeplerOrbitUtils.Epsilon && dot > -KeplerOrbitUtils.Epsilon)
				{
					return 0;
				}
				var angle = Math.Acos(dot);
				if (KeplerOrbitUtils.DotProduct(KeplerOrbitUtils.CrossProduct(ascNodeDir, SemiMajorAxisBasis), OrbitNormal) < 0)
				{
					angle = KeplerOrbitUtils.PI_2 - angle;
				}
				return angle;
			}
		}

		/// <summary>
		/// Is orbit state valid and error-free.
		/// </summary>
		/// <value>
		///   <c>true</c> if this instance is valid orbit; otherwise, <c>false</c>.
		/// </value>
		public bool IsValidOrbit
		{
			get
			{
				return Eccentricity >= 0 
					&& Period > KeplerOrbitUtils.Epsilon 
					&& AttractorDistance > KeplerOrbitUtils.Epsilon 
					&& AttractorMass > KeplerOrbitUtils.Epsilon;
			}
		}

		/// <summary>
		/// Create new orbit state without initialization. Manual orbit initialization is required.
		/// </summary>
		/// <remarks>
		/// To manually initialize orbit, fill known orbital elements and then call CalculateOrbitStateFrom... method.
		/// </remarks>
		public KeplerOrbitData()
		{
		}

		/// <summary>
		/// Create and initialize new orbit state.
		/// </summary>
		/// <param name="position">Body local position, relative to attractor.</param>
		/// <param name="velocity">Body local velocity.</param>
		/// <param name="attractorMass">Attractor mass.</param>
		/// <param name="gConst">Gravitational Constant.</param>
		public KeplerOrbitData(Vector3d position, Vector3d velocity, double attractorMass, double gConst)
		{
			this.Position = position;
			this.Velocity = velocity;
			this.AttractorMass = attractorMass;
			this.GravConst = gConst;
			CalculateOrbitStateFromOrbitalVectors();
		}

		/// <summary>
		/// Create and initialize new orbit state from orbital elements.
		/// </summary>
		/// <param name="eccentricity">Eccentricity.</param>
		/// <param name="semiMajorAxis">Main axis semi width.</param>
		/// <param name="meanAnomalyDeg">Mean anomaly in degrees.</param>
		/// <param name="inclinationDeg">Orbit inclination in degrees.</param>
		/// <param name="argOfPerifocus">Orbit argument of perifocus in degrees.</param>
		/// <param name="ascendingNodeDeg">Longitude of ascending node in degrees.</param>
		/// <param name="attractorMass">Attractor mass.</param>
		/// <param name="gConst">Gravitational constant.</param>
		public KeplerOrbitData(double eccentricity, double semiMajorAxis, double meanAnomalyDeg, double inclinationDeg, double argOfPerifocus, double ascendingNodeDeg, double attractorMass, double gConst)
		{
			this.Eccentricity = eccentricity;
			this.SemiMajorAxis = semiMajorAxis;
			if (eccentricity < 1.0)
			{
				this.SemiMinorAxis = SemiMajorAxis * Math.Sqrt(1 - Eccentricity * Eccentricity);
			}
			else
			{
				this.SemiMinorAxis = SemiMajorAxis * Math.Sqrt(Eccentricity * Eccentricity - 1);
			}
			
			var normal = EclipticNormal.normalized;
			var ascendingNode = EclipticRight.normalized;

			ascendingNodeDeg %= 360;
			if (ascendingNodeDeg > 180) ascendingNodeDeg -= 360;
			inclinationDeg %= 360;
			if (inclinationDeg > 180) inclinationDeg -= 360;
			argOfPerifocus %= 360;
			if (argOfPerifocus > 180) argOfPerifocus -= 360;

			ascendingNode = KeplerOrbitUtils.RotateVectorByAngle(ascendingNode, ascendingNodeDeg * KeplerOrbitUtils.Deg2Rad, normal).normalized;
			normal = KeplerOrbitUtils.RotateVectorByAngle(normal, inclinationDeg * KeplerOrbitUtils.Deg2Rad, ascendingNode).normalized;
			var periapsis = ascendingNode;
			periapsis = KeplerOrbitUtils.RotateVectorByAngle(periapsis, argOfPerifocus * KeplerOrbitUtils.Deg2Rad, normal).normalized;

			this.SemiMajorAxisBasis = periapsis;
			this.SemiMinorAxisBasis = KeplerOrbitUtils.CrossProduct(periapsis, normal);

			this.MeanAnomaly = meanAnomalyDeg * KeplerOrbitUtils.Deg2Rad;
			this.EccentricAnomaly = KeplerOrbitUtils.ConvertMeanToEccentricAnomaly(this.MeanAnomaly, this.Eccentricity);
			this.TrueAnomaly = KeplerOrbitUtils.ConvertEccentricToTrueAnomaly(this.EccentricAnomaly, this.Eccentricity);
			this.AttractorMass = attractorMass;
			this.GravConst = gConst;
			CalculateOrbitStateFromOrbitalElements();
		}

		/// <summary>
		/// Create and initialize new orbit state from orbital elements and main axis vectors.
		/// </summary>
		/// <param name="eccentricity">Eccentricity.</param>
		/// <param name="semiMajorAxis">Semi major axis vector.</param>
		/// <param name="semiMinorAxis">Semi minor axis vector.</param>
		/// <param name="meanAnomalyDeg">Mean anomaly in degrees.</param>
		/// <param name="attractorMass">Attractor mass.</param>
		/// <param name="gConst">Gravitational constant.</param>
		public KeplerOrbitData(double eccentricity, Vector3d semiMajorAxis, Vector3d semiMinorAxis, double meanAnomalyDeg, double attractorMass, double gConst)
		{
			this.Eccentricity = eccentricity;
			this.SemiMajorAxisBasis = semiMajorAxis.normalized;
			this.SemiMinorAxisBasis = semiMinorAxis.normalized;
			this.SemiMajorAxis = semiMajorAxis.magnitude;
			this.SemiMinorAxis = semiMinorAxis.magnitude;

			this.MeanAnomaly = meanAnomalyDeg * KeplerOrbitUtils.Deg2Rad;
			this.EccentricAnomaly = KeplerOrbitUtils.ConvertMeanToEccentricAnomaly(this.MeanAnomaly, this.Eccentricity);
			this.TrueAnomaly = KeplerOrbitUtils.ConvertEccentricToTrueAnomaly(this.EccentricAnomaly, this.Eccentricity);
			this.AttractorMass = attractorMass;
			this.GravConst = gConst;
			CalculateOrbitStateFromOrbitalElements();
		}

		[Obsolete("Use CalculateOrbitStateFromOrbitalVectors() instead.", error: true)]
		public void CalculateNewOrbitData()
		{
			CalculateOrbitStateFromOrbitalVectors();
		}

		/// <summary>
		/// Calculates full orbit state from cartesian vectors: current body position, velocity, attractor mass, and gravConstant.
		/// </summary>
		public void CalculateOrbitStateFromOrbitalVectors()
		{
			MG = AttractorMass * GravConst;
			AttractorDistance = Position.magnitude;
			Vector3d angularMomentumVector = KeplerOrbitUtils.CrossProduct(Position, Velocity);
			OrbitNormal = angularMomentumVector.normalized;
			Vector3d eccVector;
			if (OrbitNormal.sqrMagnitude < 0.99)
			{
				// If normalized vector len is not one, then it's zero.
				OrbitNormal = KeplerOrbitUtils.CrossProduct(Position, EclipticUp).normalized;
				eccVector = new Vector3d();
			}
			else
			{
				eccVector = KeplerOrbitUtils.CrossProduct(Velocity, angularMomentumVector) / MG - Position / AttractorDistance;
			}
			OrbitNormalDotEclipticNormal = KeplerOrbitUtils.DotProduct(OrbitNormal, EclipticNormal);
			FocalParameter = angularMomentumVector.sqrMagnitude / MG;
			Eccentricity = eccVector.magnitude;
			EnergyTotal = Velocity.sqrMagnitude - 2 * MG / AttractorDistance;
			SemiMinorAxisBasis = KeplerOrbitUtils.CrossProduct(angularMomentumVector, -eccVector).normalized;
			if (SemiMinorAxisBasis.sqrMagnitude < 0.99)
			{
				SemiMinorAxisBasis = KeplerOrbitUtils.CrossProduct(OrbitNormal, Position).normalized;
			}
			SemiMajorAxisBasis = KeplerOrbitUtils.CrossProduct(OrbitNormal, SemiMinorAxisBasis).normalized;
			if (Eccentricity < 1)
			{
				OrbitCompressionRatio = 1 - Eccentricity * Eccentricity;
				SemiMajorAxis = FocalParameter / OrbitCompressionRatio;
				SemiMinorAxis = SemiMajorAxis * Math.Sqrt(OrbitCompressionRatio);
				CenterPoint = -SemiMajorAxis * eccVector;
				Period = KeplerOrbitUtils.PI_2 * Math.Sqrt(Math.Pow(SemiMajorAxis, 3) / MG);
				Apoapsis = CenterPoint - SemiMajorAxisBasis * SemiMajorAxis;
				Periapsis = CenterPoint + SemiMajorAxisBasis * SemiMajorAxis;
				PeriapsisDistance = Periapsis.magnitude;
				ApoapsisDistance = Apoapsis.magnitude;
				TrueAnomaly = Vector3d.Angle(Position, SemiMajorAxisBasis) * KeplerOrbitUtils.Deg2Rad;
				if (KeplerOrbitUtils.DotProduct(KeplerOrbitUtils.CrossProduct(Position, -SemiMajorAxisBasis), OrbitNormal) < 0)
				{
					TrueAnomaly = KeplerOrbitUtils.PI_2 - TrueAnomaly;
				}
				EccentricAnomaly = KeplerOrbitUtils.ConvertTrueToEccentricAnomaly(TrueAnomaly, Eccentricity);
				MeanAnomaly = EccentricAnomaly - Eccentricity * Math.Sin(EccentricAnomaly);
			}
			else
			{
				OrbitCompressionRatio = Eccentricity * Eccentricity - 1;
				SemiMajorAxis = FocalParameter / OrbitCompressionRatio;
				SemiMinorAxis = SemiMajorAxis * Math.Sqrt(OrbitCompressionRatio);
				CenterPoint = SemiMajorAxis * eccVector;
				Period = double.PositiveInfinity;
				Apoapsis = new Vector3d(double.PositiveInfinity, double.PositiveInfinity, double.PositiveInfinity);
				Periapsis = CenterPoint - SemiMajorAxisBasis * (SemiMajorAxis);
				PeriapsisDistance = Periapsis.magnitude;
				ApoapsisDistance = double.PositiveInfinity;
				TrueAnomaly = Vector3d.Angle(Position, eccVector) * KeplerOrbitUtils.Deg2Rad;
				if (KeplerOrbitUtils.DotProduct(KeplerOrbitUtils.CrossProduct(Position, -SemiMajorAxisBasis), OrbitNormal) < 0)
				{
					TrueAnomaly = -TrueAnomaly;
				}
				EccentricAnomaly = KeplerOrbitUtils.ConvertTrueToEccentricAnomaly(TrueAnomaly, Eccentricity);
				MeanAnomaly = Math.Sinh(EccentricAnomaly) * Eccentricity - EccentricAnomaly;
			}
		}

		/// <summary>
		/// Calculates the full state of orbit from current orbital elements: eccentricity, mean anomaly, semi major and semi minor axis.
		/// </summary>
		/// <remarks>
		/// Update orbital state using known main orbital elements and basis axis vectors.
		/// Can be used for first initialization of orbit state, in this case initial data must be filled before this method call.
		/// Required initial data: eccentricity, mean anomaly, inclination, attractor mass, grav constant, all anomalies, semi minor and semi major axis vectors and magnitudes.
		/// Note that semi minor and semi major axis must be fully precalculated from inclination and argument of periapsis or another source data;
		/// </remarks>
		public void CalculateOrbitStateFromOrbitalElements()
		{
			MG = AttractorMass * GravConst;
			OrbitNormal = -KeplerOrbitUtils.CrossProduct(SemiMajorAxisBasis, SemiMinorAxisBasis).normalized;
			OrbitNormalDotEclipticNormal = KeplerOrbitUtils.DotProduct(OrbitNormal, EclipticNormal);
			if (Eccentricity < 1.0)
			{
				OrbitCompressionRatio = 1 - Eccentricity * Eccentricity;
				CenterPoint = -SemiMajorAxisBasis * SemiMajorAxis * Eccentricity;
				Period = KeplerOrbitUtils.PI_2 * Math.Sqrt(Math.Pow(SemiMajorAxis, 3) / MG);
				Apoapsis = CenterPoint - SemiMajorAxisBasis * SemiMajorAxis;
				Periapsis = CenterPoint + SemiMajorAxisBasis * SemiMajorAxis;
				PeriapsisDistance = Periapsis.magnitude;
				ApoapsisDistance = Apoapsis.magnitude;
				// All anomalies state already preset.
			}
			else
			{
				CenterPoint = SemiMajorAxisBasis * SemiMajorAxis * Eccentricity;
				Period = double.PositiveInfinity;
				Apoapsis = new Vector3d(double.PositiveInfinity, double.PositiveInfinity, double.PositiveInfinity);
				Periapsis = CenterPoint - SemiMajorAxisBasis * (SemiMajorAxis);
				PeriapsisDistance = Periapsis.magnitude;
				ApoapsisDistance = double.PositiveInfinity;
			}
			Position = GetFocalPositionAtEccentricAnomaly(EccentricAnomaly);
			double compresion = Eccentricity < 1 ? (1 - Eccentricity * Eccentricity) : (Eccentricity * Eccentricity - 1);
			FocalParameter = SemiMajorAxis * compresion;
			Velocity = GetVelocityAtTrueAnomaly(this.TrueAnomaly);
			AttractorDistance = Position.magnitude;
			EnergyTotal = Velocity.sqrMagnitude - 2 * MG / AttractorDistance;
		}

		/// <summary>
		/// Gets the velocity vector value at eccentric anomaly.
		/// </summary>
		/// <param name="eccentricAnomaly">The eccentric anomaly.</param>
		/// <returns>Velocity vector.</returns>
		public Vector3d GetVelocityAtEccentricAnomaly(double eccentricAnomaly)
		{
			return GetVelocityAtTrueAnomaly(KeplerOrbitUtils.ConvertEccentricToTrueAnomaly(eccentricAnomaly, Eccentricity));
		}

		/// <summary>
		/// Gets the velocity value at true anomaly.
		/// </summary>
		/// <param name="trueAnomaly">The true anomaly.</param>
		/// <returns>Velocity vector.</returns>
		public Vector3d GetVelocityAtTrueAnomaly(double trueAnomaly)
		{
			if (FocalParameter < KeplerOrbitUtils.Epsilon)
			{
				return new Vector3d();
			}
			double sqrtMGdivP = Math.Sqrt(AttractorMass * GravConst / FocalParameter);
			double vX = sqrtMGdivP * (Eccentricity + Math.Cos(trueAnomaly));
			double vY = sqrtMGdivP * Math.Sin(trueAnomaly);
			return -SemiMinorAxisBasis * vX - SemiMajorAxisBasis * vY;
		}

		/// <summary>
		/// Gets the central position at true anomaly.
		/// </summary>
		/// <param name="trueAnomaly">The true anomaly.</param>
		/// <returns>Position relative to orbit center.</returns>
		/// <remarks>
		/// Note: central position is not same as focal position.
		/// </remarks>
		public Vector3d GetCentralPositionAtTrueAnomaly(double trueAnomaly)
		{
			double ecc = KeplerOrbitUtils.ConvertTrueToEccentricAnomaly(trueAnomaly, Eccentricity);
			return GetCentralPositionAtEccentricAnomaly(ecc);
		}

		/// <summary>
		/// Gets the central position at eccentric anomaly.
		/// </summary>
		/// <param name="eccentricAnomaly">The eccentric anomaly.</param>
		/// <returns>Position relative to orbit center.</returns>
		/// <remarks>
		/// Note: central position is not same as focal position.
		/// </remarks>
		public Vector3d GetCentralPositionAtEccentricAnomaly(double eccentricAnomaly)
		{
			Vector3d result = Eccentricity < 1 ?
				new Vector3d(Math.Sin(eccentricAnomaly) * SemiMinorAxis, -Math.Cos(eccentricAnomaly) * SemiMajorAxis) :
				new Vector3d(Math.Sinh(eccentricAnomaly) * SemiMinorAxis, Math.Cosh(eccentricAnomaly) * SemiMajorAxis);
			return -SemiMinorAxisBasis * result.x - SemiMajorAxisBasis * result.y;
		}

		/// <summary>
		/// Gets the focal position at eccentric anomaly.
		/// </summary>
		/// <param name="eccentricAnomaly">The eccentric anomaly.</param>
		/// <returns>Position relative to attractor (focus).</returns>
		public Vector3d GetFocalPositionAtEccentricAnomaly(double eccentricAnomaly)
		{
			return GetCentralPositionAtEccentricAnomaly(eccentricAnomaly) + CenterPoint;
		}

		/// <summary>
		/// Gets the focal position at true anomaly.
		/// </summary>
		/// <param name="trueAnomaly">The true anomaly.</param>
		/// <returns>Position relative to attractor (focus).</returns>
		public Vector3d GetFocalPositionAtTrueAnomaly(double trueAnomaly)
		{
			return GetCentralPositionAtTrueAnomaly(trueAnomaly) + CenterPoint;
		}

		/// <summary>
		/// Gets the central position.
		/// </summary>
		/// <returns>Position relative to orbit center.</returns>
		/// <remarks>
		/// Note: central position is not same as focal position.
		/// </remarks>
		public Vector3d GetCentralPosition()
		{
			return Position - CenterPoint;
		}

		/// <summary>
		/// Gets calculated orbit points with defined precision.
		/// </summary>
		/// <param name="pointsCount">The points count.</param>
		/// <param name="maxDistance">The maximum distance of points.</param>
		/// <returns>Array of orbit curve points.</returns>
		public Vector3d[] GetOrbitPoints(int pointsCount = 50, double maxDistance = 1000d)
		{
			return GetOrbitPoints(pointsCount, new Vector3d(), maxDistance);
		}

		/// <summary>
		/// Gets calculated orbit points with defined precision.
		/// </summary>
		/// <param name="pointsCount">The points count.</param>
		/// <param name="origin">The origin.</param>
		/// <param name="maxDistance">The maximum distance.</param>
		/// <returns>Array of orbit curve points.</returns>
		public Vector3d[] GetOrbitPoints(int pointsCount, Vector3d origin, double maxDistance = 1000d)
		{
			if (pointsCount < 2)
			{
				return new Vector3d[0];
			}
			Vector3d[] result = new Vector3d[pointsCount];
			if (Eccentricity < 1)
			{
				if (ApoapsisDistance < maxDistance)
				{
					for (int i = 0; i < pointsCount; i++)
					{
						result[i] = GetFocalPositionAtEccentricAnomaly(i * KeplerOrbitUtils.PI_2 / (pointsCount - 1d)) + origin;
					}
				}
				else
				{
					double maxAngle = KeplerOrbitUtils.CalcTrueAnomalyForDistance(maxDistance, Eccentricity, SemiMajorAxis);
					for (int i = 0; i < pointsCount; i++)
					{
						result[i] = GetFocalPositionAtTrueAnomaly(-maxAngle + i * 2d * maxAngle / (pointsCount - 1)) + origin;
					}
				}
			}
			else
			{
				if (maxDistance < PeriapsisDistance)
				{
					return new Vector3d[0];
				}
				double maxAngle = KeplerOrbitUtils.CalcTrueAnomalyForDistance(maxDistance, Eccentricity, SemiMajorAxis);

				for (int i = 0; i < pointsCount; i++)
				{
					result[i] = GetFocalPositionAtTrueAnomaly(-maxAngle + i * 2d * maxAngle / (pointsCount - 1)) + origin;
				}
			}
			return result;
		}

		/// <summary>
		/// Gets calculated orbit points with defined precision.
		/// </summary>
		/// <param name="pointsCount">The points count.</param>
		/// <param name="maxDistance">The maximum distance.</param>
		/// <returns>Array of orbit curve points.</returns>
		public Vector3[] GetOrbitPoints(int pointsCount = 50, float maxDistance = 1000f)
		{
			return GetOrbitPoints(pointsCount, new Vector3(), maxDistance);
		}

		/// <summary>
		/// Gets calculated orbit points with defined precision.
		/// </summary>
		/// <param name="pointsCount">The points count.</param>
		/// <param name="origin">The origin.</param>
		/// <param name="maxDistance">The maximum distance.</param>
		/// <returns>Array of orbit curve points.</returns>
		public Vector3[] GetOrbitPoints(int pointsCount, Vector3 origin, float maxDistance = 1000f)
		{
			if (pointsCount < 2)
			{
				return new Vector3[0];
			}
			Vector3[] result = new Vector3[pointsCount];
			if (Eccentricity < 1)
			{
				if (ApoapsisDistance < maxDistance)
				{
					for (int i = 0; i < pointsCount; i++)
					{
						result[i] = (Vector3)GetFocalPositionAtEccentricAnomaly(i * KeplerOrbitUtils.PI_2 / (pointsCount - 1d)) + origin;
					}
				}
				else
				{
					double maxAngle = KeplerOrbitUtils.CalcTrueAnomalyForDistance(maxDistance, Eccentricity, SemiMajorAxis);
					for (int i = 0; i < pointsCount; i++)
					{
						result[i] = (Vector3)GetFocalPositionAtTrueAnomaly(-maxAngle + i * 2d * maxAngle / (pointsCount - 1)) + origin;
					}
				}
			}
			else
			{
				if (maxDistance < PeriapsisDistance)
				{
					return new Vector3[0];
				}
				double maxAngle = KeplerOrbitUtils.CalcTrueAnomalyForDistance(maxDistance, Eccentricity, SemiMajorAxis);

				for (int i = 0; i < pointsCount; i++)
				{
					result[i] = (Vector3)GetFocalPositionAtTrueAnomaly(-maxAngle + i * 2d * maxAngle / (pointsCount - 1)) + origin;
				}
			}
			return result;
		}

		/// <summary>
		/// Gets the orbit points without unnecessary memory alloc for resulting array.
		/// However, memory allocation may occur if resulting array has not correct lenght.
		/// </summary>
		/// <param name="orbitPoints">The orbit points.</param>
		/// <param name="pointsCount">The points count.</param>
		/// <param name="origin">The origin.</param>
		/// <param name="maxDistance">The maximum distance.</param>
		public void GetOrbitPointsNoAlloc(ref Vector3[] orbitPoints, int pointsCount, Vector3 origin, float maxDistance = 1000f)
		{
			if (pointsCount < 2)
			{
				orbitPoints = new Vector3[0];
				return;
			}
			if (Eccentricity < 1)
			{
				if (orbitPoints == null || orbitPoints.Length != pointsCount)
				{
					orbitPoints = new Vector3[pointsCount];
				}
				if (ApoapsisDistance < maxDistance)
				{
					for (int i = 0; i < pointsCount; i++)
					{
						orbitPoints[i] = (Vector3)GetFocalPositionAtEccentricAnomaly(i * KeplerOrbitUtils.PI_2 / (pointsCount - 1d)) + origin;
					}
				}
				else
				{
					double maxAngle = KeplerOrbitUtils.CalcTrueAnomalyForDistance(maxDistance, Eccentricity, SemiMajorAxis);
					for (int i = 0; i < pointsCount; i++)
					{
						orbitPoints[i] = (Vector3)GetFocalPositionAtTrueAnomaly(-maxAngle + i * 2d * maxAngle / (pointsCount - 1)) + origin;
					}
				}
			}
			else
			{
				if (maxDistance < PeriapsisDistance)
				{
					orbitPoints = new Vector3[0];
					return;
				}
				if (orbitPoints == null || orbitPoints.Length != pointsCount)
				{
					orbitPoints = new Vector3[pointsCount];
				}
				double maxAngle = KeplerOrbitUtils.CalcTrueAnomalyForDistance(maxDistance, Eccentricity, SemiMajorAxis);

				for (int i = 0; i < pointsCount; i++)
				{
					orbitPoints[i] = (Vector3)GetFocalPositionAtTrueAnomaly(-maxAngle + i * 2d * maxAngle / (pointsCount - 1)) + origin;
				}
			}
		}

		/// <summary>
		/// Gets the ascending node of orbit.
		/// </summary>
		/// <param name="asc">The asc.</param>
		/// <returns><c>true</c> if ascending node exists, otherwise <c>false</c></returns>
		public bool GetAscendingNode(out Vector3 asc)
		{
			Vector3d v;
			if (GetAscendingNode(out v))
			{
				asc = (Vector3)v;
				return true;
			}
			asc = new Vector3();
			return false;
		}

		/// <summary>
		/// Gets the ascending node of orbit.
		/// </summary>
		/// <param name="asc">The asc.</param>
		/// <returns><c>true</c> if ascending node exists, otherwise <c>false</c></returns>
		public bool GetAscendingNode(out Vector3d asc)
		{
			Vector3d ascNodeDir = KeplerOrbitUtils.CrossProduct(OrbitNormal, EclipticNormal);
			bool s = KeplerOrbitUtils.DotProduct(KeplerOrbitUtils.CrossProduct(ascNodeDir, SemiMajorAxisBasis), OrbitNormal) >= 0;
			double ecc = 0d;
			double trueAnom = Vector3d.Angle(ascNodeDir, CenterPoint) * KeplerOrbitUtils.Deg2Rad;
			if (Eccentricity < 1)
			{
				double cosT = Math.Cos(trueAnom);
				ecc = Math.Acos((Eccentricity + cosT) / (1d + Eccentricity * cosT));
				if (!s)
				{
					ecc = KeplerOrbitUtils.PI_2 - ecc;
				}
			}
			else
			{
				trueAnom = Vector3d.Angle(-ascNodeDir, CenterPoint) * KeplerOrbitUtils.Deg2Rad;
				if (trueAnom >= Math.Acos(-1d / Eccentricity))
				{
					asc = new Vector3d();
					return false;
				}
				double cosT = Math.Cos(trueAnom);
				ecc = KeplerOrbitUtils.Acosh((Eccentricity + cosT) / (1 + Eccentricity * cosT)) * (!s ? -1 : 1);
			}
			asc = GetFocalPositionAtEccentricAnomaly(ecc);
			return true;
		}

		/// <summary>
		/// Gets the descending node of orbit.
		/// </summary>
		/// <param name="desc">The desc.</param>
		/// <returns><c>true</c> if descending node exists, otherwise <c>false</c></returns>
		public bool GetDescendingNode(out Vector3 desc)
		{
			Vector3d v;
			if (GetDescendingNode(out v))
			{
				desc = (Vector3)v;
				return true;
			}
			desc = new Vector3();
			return false;
		}

		/// <summary>
		/// Gets the descending node orbit.
		/// </summary>
		/// <param name="desc">The desc.</param>
		/// <returns><c>true</c> if descending node exists, otherwise <c>false</c></returns>
		public bool GetDescendingNode(out Vector3d desc)
		{
			Vector3d norm = KeplerOrbitUtils.CrossProduct(OrbitNormal, EclipticNormal);
			bool s = KeplerOrbitUtils.DotProduct(KeplerOrbitUtils.CrossProduct(norm, SemiMajorAxisBasis), OrbitNormal) < 0;
			double ecc = 0d;
			double trueAnom = Vector3d.Angle(norm, -CenterPoint) * KeplerOrbitUtils.Deg2Rad;
			if (Eccentricity < 1)
			{
				double cosT = Math.Cos(trueAnom);
				ecc = Math.Acos((Eccentricity + cosT) / (1d + Eccentricity * cosT));
				if (s)
				{
					ecc = KeplerOrbitUtils.PI_2 - ecc;
				}
			}
			else
			{
				trueAnom = Vector3d.Angle(norm, CenterPoint) * KeplerOrbitUtils.Deg2Rad;
				if (trueAnom >= Math.Acos(-1d / Eccentricity))
				{
					desc = new Vector3d();
					return false;
				}
				double cosT = Math.Cos(trueAnom);
				ecc = KeplerOrbitUtils.Acosh((Eccentricity + cosT) / (1 + Eccentricity * cosT)) * (s ? -1 : 1);
			}
			desc = GetFocalPositionAtEccentricAnomaly(ecc);
			return true;
		}

		/// <summary>
		/// Updates the kepler orbit state by defined deltatime.
		/// Orbit main parameters will remains unchanged, but all anomalies will progress in time.
		/// </summary>
		/// <param name="deltaTime">The delta time.</param>
		public void UpdateOrbitDataByTime(double deltaTime)
		{
			UpdateOrbitAnomaliesByTime(deltaTime);
			SetPositionByCurrentAnomaly();
			SetVelocityByCurrentAnomaly();
		}

		/// <summary>
		/// Updates the value of orbital anomalies by defined deltatime.
		/// </summary>
		/// <param name="deltaTime">The delta time.</param>
		/// <remarks>
		/// Only anomalies values will be changed. 
		/// Position and velocity states needs to be updated too after this method call.
		/// </remarks>
		public void UpdateOrbitAnomaliesByTime(double deltaTime)
		{
			if (Eccentricity < 1)
			{
				if (Period > KeplerOrbitUtils.Epsilon)
				{
					MeanAnomaly += KeplerOrbitUtils.PI_2 * deltaTime / Period;
				}
				MeanAnomaly %= KeplerOrbitUtils.PI_2;
				if (MeanAnomaly < 0)
				{
					MeanAnomaly = KeplerOrbitUtils.PI_2 - MeanAnomaly;
				}
				EccentricAnomaly = KeplerOrbitUtils.KeplerSolver(MeanAnomaly, Eccentricity);
				double cosE = Math.Cos(EccentricAnomaly);
				TrueAnomaly = Math.Acos((cosE - Eccentricity) / (1 - Eccentricity * cosE));
				if (MeanAnomaly > Math.PI)
				{
					TrueAnomaly = KeplerOrbitUtils.PI_2 - TrueAnomaly;
				}
				if (double.IsNaN(MeanAnomaly) || double.IsInfinity(MeanAnomaly))
				{
					Debug.Log("KeplerOrbitData: NaN(INF) MEAN ANOMALY"); //little paranoya
					Debug.Break();
				}
				if (double.IsNaN(EccentricAnomaly) || double.IsInfinity(EccentricAnomaly))
				{
					Debug.Log("KeplerOrbitData: NaN(INF) ECC ANOMALY");
					Debug.Break();
				}
				if (double.IsNaN(TrueAnomaly) || double.IsInfinity(TrueAnomaly))
				{
					Debug.Log("KeplerOrbitData: NaN(INF) TRUE ANOMALY");
					Debug.Break();
				}
			}
			else
			{
				double n = Math.Sqrt(AttractorMass * GravConst / Math.Pow(SemiMajorAxis, 3));
				MeanAnomaly = MeanAnomaly + n * deltaTime;
				EccentricAnomaly = KeplerOrbitUtils.KeplerSolverHyperbolicCase(MeanAnomaly, Eccentricity);
				TrueAnomaly = Math.Atan2(Math.Sqrt(Eccentricity * Eccentricity - 1.0) * Math.Sinh(EccentricAnomaly), Eccentricity - Math.Cosh(EccentricAnomaly));
			}
		}

		/// <summary>
		/// Updates position from anomaly state.
		/// </summary>
		public void SetPositionByCurrentAnomaly()
		{
			Position = GetFocalPositionAtEccentricAnomaly(EccentricAnomaly);
		}

		/// <summary>
		/// Sets orbit velocity, calculated by current anomaly.
		/// </summary>
		public void SetVelocityByCurrentAnomaly()
		{
			Velocity = GetVelocityAtEccentricAnomaly(EccentricAnomaly);
		}

		/// <summary>
		/// Sets the eccentricity and updates all corresponding orbit state values.
		/// </summary>
		/// <param name="e">The new eccentricity value.</param>
		public void SetEccentricity(double e)
		{
			if (!IsValidOrbit)
			{
				return;
			}
			e = Math.Abs(e);
			double periapsis = PeriapsisDistance; // Periapsis remains constant
			Eccentricity = e;
			double compresion = Eccentricity < 1 ? (1 - Eccentricity * Eccentricity) : (Eccentricity * Eccentricity - 1);
			SemiMajorAxis = Math.Abs(periapsis / (1 - Eccentricity));
			FocalParameter = SemiMajorAxis * compresion;
			SemiMinorAxis = SemiMajorAxis * Math.Sqrt(compresion);
			CenterPoint = SemiMajorAxis * Math.Abs(Eccentricity) * SemiMajorAxisBasis;
			if (Eccentricity < 1)
			{
				EccentricAnomaly = KeplerOrbitUtils.KeplerSolver(MeanAnomaly, Eccentricity);
				double cosE = Math.Cos(EccentricAnomaly);
				TrueAnomaly = Math.Acos((cosE - Eccentricity) / (1 - Eccentricity * cosE));
				if (MeanAnomaly > Math.PI)
				{
					TrueAnomaly = KeplerOrbitUtils.PI_2 - TrueAnomaly;
				}
			}
			else
			{
				EccentricAnomaly = KeplerOrbitUtils.KeplerSolverHyperbolicCase(MeanAnomaly, Eccentricity);
				TrueAnomaly = Math.Atan2(Math.Sqrt(Eccentricity * Eccentricity - 1) * Math.Sinh(EccentricAnomaly), Eccentricity - Math.Cosh(EccentricAnomaly));
			}
			SetVelocityByCurrentAnomaly();
			SetPositionByCurrentAnomaly();

			CalculateOrbitStateFromOrbitalVectors();
		}

		/// <summary>
		/// Sets the mean anomaly and updates all other anomalies.
		/// </summary>
		/// <param name="m">The m.</param>
		public void SetMeanAnomaly(double m)
		{
			if (!IsValidOrbit)
			{
				return;
			}
			MeanAnomaly = m % KeplerOrbitUtils.PI_2;
			if (Eccentricity < 1)
			{
				if (MeanAnomaly < 0)
				{
					MeanAnomaly += KeplerOrbitUtils.PI_2;
				}
				EccentricAnomaly = KeplerOrbitUtils.KeplerSolver(MeanAnomaly, Eccentricity);
				TrueAnomaly = KeplerOrbitUtils.ConvertEccentricToTrueAnomaly(EccentricAnomaly, Eccentricity);
			}
			else
			{
				EccentricAnomaly = KeplerOrbitUtils.KeplerSolverHyperbolicCase(MeanAnomaly, Eccentricity);
				TrueAnomaly = KeplerOrbitUtils.ConvertEccentricToTrueAnomaly(EccentricAnomaly, Eccentricity);
			}
			SetPositionByCurrentAnomaly();
			SetVelocityByCurrentAnomaly();
		}

		/// <summary>
		/// Sets the true anomaly and updates all other anomalies.
		/// </summary>
		/// <param name="t">The t.</param>
		public void SetTrueAnomaly(double t)
		{
			if (!IsValidOrbit)
			{
				return;
			}
			t %= KeplerOrbitUtils.PI_2;

			if (Eccentricity < 1)
			{
				if (t < 0)
				{
					t += KeplerOrbitUtils.PI_2;
				}
				EccentricAnomaly = KeplerOrbitUtils.ConvertTrueToEccentricAnomaly(t, Eccentricity);
				MeanAnomaly = EccentricAnomaly - Eccentricity * Math.Sin(EccentricAnomaly);
			}
			else
			{
				EccentricAnomaly = KeplerOrbitUtils.ConvertTrueToEccentricAnomaly(t, Eccentricity);
				MeanAnomaly = Math.Sinh(EccentricAnomaly) * Eccentricity - EccentricAnomaly;
			}
			SetPositionByCurrentAnomaly();
			SetVelocityByCurrentAnomaly();
		}

		/// <summary>
		/// Sets the eccentric anomaly and updates all other anomalies.
		/// </summary>
		/// <param name="e">The e.</param>
		public void SetEccentricAnomaly(double e)
		{
			if (!IsValidOrbit)
			{
				return;
			}
			e %= KeplerOrbitUtils.PI_2;
			EccentricAnomaly = e;
			if (Eccentricity < 1)
			{
				if (e < 0)
				{
					e = KeplerOrbitUtils.PI_2 + e;
				}
				EccentricAnomaly = e;
				TrueAnomaly = KeplerOrbitUtils.ConvertEccentricToTrueAnomaly(e, Eccentricity);
				MeanAnomaly = EccentricAnomaly - Eccentricity * Math.Sin(EccentricAnomaly);
			}
			else
			{
				TrueAnomaly = KeplerOrbitUtils.ConvertEccentricToTrueAnomaly(e, Eccentricity);
				MeanAnomaly = Math.Sinh(EccentricAnomaly) * Eccentricity - EccentricAnomaly;
			}
			SetPositionByCurrentAnomaly();
			SetVelocityByCurrentAnomaly();
		}

		/// <summary>
		/// Rotates the relative position and velocity by same quaternion.
		/// </summary>
		/// <param name="rotation">The rotation.</param>
		public void RotateOrbit(Quaternion rotation)
		{
			Position = new Vector3d(rotation * ((Vector3)Position));
			Velocity = new Vector3d(rotation * ((Vector3)Velocity));
			CalculateOrbitStateFromOrbitalVectors();
		}

		public object Clone()
		{
			return MemberwiseClone();
		}

		public KeplerOrbitData CloneOrbit()
		{
			return (KeplerOrbitData)MemberwiseClone();
		}
	}
}
