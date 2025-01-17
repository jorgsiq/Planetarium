﻿using System;
using System.Linq;
using System.Collections.Generic;
using UnityEngine;

namespace SimpleKeplerOrbits.Examples
{
	/// <summary>
	/// Controller for spawning bodies by orbital elements values. Supports JPL database input.
	/// </summary>
	public class JPLOrbitsLoader : MonoBehaviour
	{
		[Serializable]
		public class JPLListContainer
		{
			public JPLElementsData[] OrbitsData = new JPLElementsData[0];
		}

		/// <summary>
		/// Data container for single body orbit.
		/// </summary>
		[Serializable]
		public class JPLElementsData
		{
			public string BodyName;
			public string AttractorName;
			public float AttractorMass;
			/// <summary>
			/// Eccentricity.
			/// </summary>
			[Tooltip("Eccentricity")]
			public double EC;

			/// <summary>
			/// Inclination (degrees).
			/// </summary>
			[Tooltip("Inclination")]
			public double IN;

			/// <summary>
			/// Longitude of Ascending Node (degrees).
			/// </summary>
			[Tooltip("Ascending node longitude")]
			public double OM;

			/// <summary>
			/// Argument of perifocus (degrees).
			/// </summary>
			[Tooltip("Argument of periapsis")]
			public double W;

			/// <summary>
			/// Mean anomaly (degrees).
			/// </summary>
			[Tooltip("Mean anomaly")]
			public double MA;

			/// <summary>
			/// Semi-major axis (au).
			/// </summary>
			[Tooltip("Semi-major axis")]
			public double A;
		}

		public enum LoadingType
		{
			Json,
			Scene,
		}

		public LoadingType LoadingDataSource;

		/// <summary>
		/// Gravitational constant. In this context plays role of speed muliplier.
		/// </summary>
		public double GConstant = 100;

		/// <summary>
		/// Scale multiplier: world units per 1 au.
		/// </summary>
		public float UnitsPerAU = 1f;

		public KeplerOrbitMover BodyTemplate;

		public TextAsset JsonData;

		public JPLElementsData[] ElementsTable;

		private void Start()
		{
			switch (LoadingDataSource)
			{
				case LoadingType.Json:
					SpawnFromJsonData();
					break;
				case LoadingType.Scene:
					SpawnFromSceneData();
					break;
			}
		}

		[ContextMenu("Spawn from json")]
		private void SpawnFromJsonData()
		{
			if (JsonData != null)
			{
				var data = JsonUtility.FromJson<JPLListContainer>(JsonData.text);
				if (data != null)
				{
					SpawnAll(data.OrbitsData);
				}
			}
		}

		[ContextMenu("Spawn from component")]
		private void SpawnFromSceneData()
		{
			SpawnAll(ElementsTable);
		}

		[ContextMenu("Destroy all bodies")]
		private void DestroyAll()
		{
			var instances = GameObject.FindObjectsOfType<KeplerOrbitMover>();
			for (int i = 0; i < instances.Length; i++)
			{
				if (instances[i].gameObject.activeSelf)
				{
					if (Application.isPlaying)
					{
						Destroy(instances[i].gameObject);
					}
					else
					{
						DestroyImmediate(instances[i].gameObject);
					}
				}
			}
		}
		
		private void SpawnAll(JPLElementsData[] inputData)
		{
			// Spawn all bodies in multiple passes to resolve parent-child connections.
			// All spawned bodies are instantiated from signle template, which has no visual components attached,
			// because this example is designed only for simplest orbits loading process demonstration.

			if (inputData == null || inputData.Length == 0) return;
			List<JPLElementsData> spawnOrder = new List<JPLElementsData>(inputData);
			List<KeplerOrbitMover> spawnedInstances = new List<KeplerOrbitMover>();

			if (GameObject.FindObjectsOfType<Transform>().Length > 5)
			{
				Debug.Log("Warning! too many object on scene. Don't forget to remove old spawned object before spawning new pass");
			}

			bool isAnySpawned = true;
			while (spawnOrder.Count > 0 && isAnySpawned)
			{
				isAnySpawned = false;
				for (int i = 0; i < spawnOrder.Count; i++)
				{
					var attractorName = spawnOrder[0].AttractorName != null ? spawnOrder[0].AttractorName.Trim() : "";
					bool isAttractorSpawned = 
						string.IsNullOrEmpty(attractorName)
						?true
						:spawnedInstances.Any(s => s.name == attractorName);
					if (isAttractorSpawned)
					{
						KeplerOrbitMover attractor = string.IsNullOrEmpty(attractorName)
							? null
							: spawnedInstances.First(s => s.name == attractorName);
						KeplerOrbitMover body = Instantiate(BodyTemplate, parent: attractor == null ? null : attractor.transform);
						if (!string.IsNullOrEmpty(spawnOrder[0].BodyName))
						{
							body.name = spawnOrder[0].BodyName.Trim();
						}
						if (attractor != null)
						{
							body.AttractorSettings.AttractorMass = spawnOrder[0].AttractorMass;
						}
						body.AttractorSettings.GravityConstant = (float)GConstant;
						body.AttractorSettings.AttractorObject = attractor == null ? null : attractor.transform;
						body.OrbitData = new KeplerOrbitData(
							eccentricity: spawnOrder[i].EC,
							semiMajorAxis: spawnOrder[i].A * UnitsPerAU,
							meanAnomalyDeg: spawnOrder[i].MA,
							inclinationDeg: spawnOrder[i].IN,
							argOfPerifocus: spawnOrder[i].W,
							ascendingNodeDeg: spawnOrder[i].OM,
							attractorMass: body.AttractorSettings.AttractorMass,
							gConst: GConstant);
						if (attractor != null)
						{
							body.ForceUpdateViewFromInternalState();
						}
						else
						{
							body.enabled = false;
						}
						body.gameObject.SetActive(true);
						spawnOrder.RemoveAt(0);
						i--;
						spawnedInstances.Add(body);
						isAnySpawned = true;
					}
					else
					{
						// If attractor not spawned yet, then wait for next spawn cycle pass.
					}
				}
			}
			if (!isAnySpawned && spawnOrder.Count > 0)
			{
				Debug.LogError("Couldn't spawn " + spawnOrder.Count + " because assigned attractor was not found");
			}
		}
	}
}