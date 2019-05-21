using System;
using System.Linq;
using System.Collections.Generic;
using UnityEngine;

namespace SimpleKeplerOrbits.Examples
{
	public class JPLOrbitsLoader : MonoBehaviour
	{
		[Serializable]
		public class JPLListContainer
		{
			public JPLElementsData[] OrbitsData = new JPLElementsData[0];
		}

		//data container por órbita única
		[Serializable]
		public class JPLElementsData
		{
			public string BodyName;
			public string AttractorName;
			public float AttractorMass;

			//ecentricidade
			[Tooltip("Eccentricity")]
			public double EC;


			//inclinação
			[Tooltip("Inclination")]
			public double IN;

			[Tooltip("Ascending node longitude")]
			public double OM;

			
			[Tooltip("Argument of periapsis")]
			public double W;

			[Tooltip("Mean anomaly")]
			public double MA;

			[Tooltip("Semi-major axis")]
			public double A;
		}

		public enum LoadingType
		{
			Json,
			Scene,
		}

		public LoadingType LoadingDataSource;

        //constante gravitacional
		public double GConstant = 100;

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