using System.Collections;
using UnityEngine;

namespace SimpleKeplerOrbits
{

    //componente para mover objeto no caminho elíptico ou hiperbólico ao redor do corpo do sol
	[ExecuteInEditMode]
	public class KeplerOrbitMover : MonoBehaviour
	{

        //a referência de objeto(sol) deve ser atribuída ou o movimento da órbita não funcionará
		public AttractorData AttractorSettings = new AttractorData();

        //atribuir objeto e usá-lo como controle de velocidade na exibição da cena
        [Tooltip("The velocity handle object. Assign object and use it as velocity control handle in scene view.")]
		public Transform VelocityHandle;


		//velocity handle parâmetro do tamanho de escala
		[Range(0f, 10f)]
		[Tooltip("Velocity handle scale parameter.")]
		public float VelocityHandleLenghtScale = 0f;

		// multiplicador da escala de tempo
		[Tooltip("The time scale multiplier.")]
		public float TimeScale = 1f;

		[Header("Orbit state details:")]
		[Tooltip("Internal state of orbit.")]
		public KeplerOrbitData OrbitData = new KeplerOrbitData();

		public bool LockOrbitEditing = false;

#if UNITY_EDITOR
	
		private bool _debugErrorDisplayed = false;
#endif

		private bool IsReferencesAsigned
		{
			get
			{
				return AttractorSettings != null && AttractorSettings.AttractorObject != null;
			}
		}

		private void OnEnable()
		{
			ForceUpdateOrbitData();
#if UNITY_EDITOR
			if (!Application.isPlaying)
			{
				return;
			}
#endif
			StartCoroutine(OrbitUpdateLoop());
		}

		private void Update()
		{
			if (IsReferencesAsigned)
			{
				if (!LockOrbitEditing)
				{
					Vector3d position = new Vector3d(transform.position - AttractorSettings.AttractorObject.position);

					bool velocityHandleChanged = false;
					if (VelocityHandle != null)
					{
						Vector3 velocity = GetVelocityHandleDisplayedVelocity();
						if (velocity != (Vector3)OrbitData.Velocity)
						{
							velocityHandleChanged = true;
						}
					}
					if ((Vector3)position != (Vector3)OrbitData.Position ||
						velocityHandleChanged ||
						OrbitData.GravConst != AttractorSettings.GravityConstant ||
						OrbitData.AttractorMass != AttractorSettings.AttractorMass)
					{
						ForceUpdateOrbitData();
					}
				}
			}
			else
			{
#if UNITY_EDITOR
				if (AttractorSettings.AttractorObject == null)
				{
					if (!_debugErrorDisplayed)
					{
						_debugErrorDisplayed = true;
						if (Application.isPlaying)
						{
							Debug.LogError("KeplerMover: Attractor reference not asigned", context: gameObject);
						}
						else
						{
							Debug.Log("KeplerMover: Attractor reference not asigned", context: gameObject);
						}
					}
				}
				else
				{
					_debugErrorDisplayed = false;
				}
#endif
			}
		}



		private IEnumerator OrbitUpdateLoop()
		{
			while (true)
			{
				if (IsReferencesAsigned)
				{
					if (!OrbitData.IsValidOrbit)
					{
						//corrige a órbita
						OrbitData.CalculateOrbitStateFromOrbitalVectors();
					}

					if (OrbitData.IsValidOrbit)
					{
						OrbitData.UpdateOrbitDataByTime(Time.deltaTime * TimeScale);
						ForceUpdateViewFromInternalState();
					}
				}
				yield return null;
			}
		}

		public void CreateNewOrbitFromPositionAndVelocity(Vector3 relativePosition, Vector3 velocity)
		{
			if (IsReferencesAsigned)
			{
				OrbitData.Position = new Vector3d(relativePosition);
				OrbitData.Velocity = new Vector3d(velocity);
				OrbitData.CalculateOrbitStateFromOrbitalVectors();
				ForceUpdateViewFromInternalState();
			}
		}


		[ContextMenu("Update transform from orbit state")]
		public void ForceUpdateViewFromInternalState()
		{
			transform.position = AttractorSettings.AttractorObject.position + (Vector3)OrbitData.Position;
			ForceUpdateVelocityHandleFromInternalState();
		}

		public void ForceUpdateVelocityHandleFromInternalState()
		{
			if (VelocityHandle != null)
			{
				Vector3 velocityRelativePosition = (Vector3)OrbitData.Velocity;
				if (VelocityHandleLenghtScale > 0 && !float.IsNaN(VelocityHandleLenghtScale) && !float.IsInfinity(VelocityHandleLenghtScale))
				{
					velocityRelativePosition *= VelocityHandleLenghtScale;
				}
				VelocityHandle.position = transform.position + velocityRelativePosition;
			}
		}

		public Vector3 GetVelocityHandleDisplayedVelocity()
		{
			if (VelocityHandle != null)
			{
				Vector3 velocity = VelocityHandle.position - transform.position;
				if (VelocityHandleLenghtScale > 0 && !float.IsNaN(VelocityHandleLenghtScale) && !float.IsInfinity(VelocityHandleLenghtScale))
				{
					velocity /= VelocityHandleLenghtScale;
				}
				return velocity;
			}
			return new Vector3();
		}

		[ContextMenu("Update Orbit data from current vectors")]
		public void ForceUpdateOrbitData()
		{
			if (IsReferencesAsigned)
			{
				OrbitData.AttractorMass = AttractorSettings.AttractorMass;
				OrbitData.GravConst = AttractorSettings.GravityConstant;
				OrbitData.Position = new Vector3d(transform.position - AttractorSettings.AttractorObject.position);
				if (VelocityHandle != null)
				{
					Vector3 velocity = GetVelocityHandleDisplayedVelocity();
					OrbitData.Velocity = new Vector3d(velocity);
				}
				OrbitData.CalculateOrbitStateFromOrbitalVectors();
			}
		}


        //altera o vetor de velocidade da órbita para coincidir com a órbita circular
        [ContextMenu("Circularize orbit")]
		public void SetAutoCircleOrbit()
		{
			if (IsReferencesAsigned)
			{
				OrbitData.Velocity = KeplerOrbitUtils.CalcCircleOrbitVelocity(Vector3d.zero, OrbitData.Position, OrbitData.AttractorMass, 1f, OrbitData.OrbitNormal, OrbitData.GravConst);
				OrbitData.CalculateOrbitStateFromOrbitalVectors();
				ForceUpdateVelocityHandleFromInternalState();
			}
		}

		[ContextMenu("Inverse velocity")]
		public void InverseVelocity()
		{
			if (IsReferencesAsigned)
			{
				OrbitData.Velocity = -OrbitData.Velocity;
				OrbitData.CalculateOrbitStateFromOrbitalVectors();
				ForceUpdateVelocityHandleFromInternalState();
			}
		}

		[ContextMenu("Inverse position")]
		public void InversePositionRelativeToAttractor()
		{
			if (IsReferencesAsigned)
			{
				OrbitData.Position = -OrbitData.Position;
				OrbitData.CalculateOrbitStateFromOrbitalVectors();
				ForceUpdateVelocityHandleFromInternalState();
			}
		}

		[ContextMenu("Inverse velocity and position")]
		public void InverseOrbit()
		{
			if (IsReferencesAsigned)
			{
				OrbitData.Velocity = -OrbitData.Velocity;
				OrbitData.Position = -OrbitData.Position;
				OrbitData.CalculateOrbitStateFromOrbitalVectors();
				ForceUpdateVelocityHandleFromInternalState();
			}
		}

		[ContextMenu("Reset orbit")]
		public void ResetOrbit()
		{
			OrbitData = new KeplerOrbitData();
		}
	}
}