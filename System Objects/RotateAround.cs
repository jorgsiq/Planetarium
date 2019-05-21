using UnityEngine;
using System.Collections;

public class RotateAround : MonoBehaviour {
    
    public Transform centerMass;
    public Transform revealingLight;

    public float rotationAroundPlanetSpeed = 1.0f;
    public float rotationAroundSunDays = 1.0f;
    public float defaultEarthYear = 365.25f;
    public float satelliteOrbitDistance = 1.8f;
    public float planetSunDistance = 100.0f;
    public float planetSpeedRotation = 1.0f;

    private GlobalValues globalValuesScript;

    void Start () {
        globalValuesScript = Camera.main.GetComponent<GlobalValues>();

        rotationAroundPlanetSpeed = rotationAroundSunDays / defaultEarthYear;

        if(revealingLight != null)
        {
            revealingLight.transform.LookAt(transform.position);
        }

    }
	

	void Update () {
    

        transform.RotateAround(centerMass.position, Vector3.up, Time.deltaTime * (defaultEarthYear / rotationAroundSunDays) * (globalValuesScript.globalPlanetRotationAroundSun) * Time.deltaTime);

        transform.Rotate(-Vector3.up * Time.deltaTime * planetSpeedRotation * globalValuesScript.globalPlanetRotationAroundSun);
    }
    
}
