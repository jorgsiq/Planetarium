using UnityEngine;
using System.Collections;
using UnityEngine.UI;

public class PlanetSwitch : MonoBehaviour {
	
    GameObject MainCamera;
    
    const float optimalDistance = 2.3f;

	public GUISkin guiSkin;

	public string planetName;

	public float baseDistance = 2.5f;
	public float baseScrollSpeed = 0.1f;

	public Light SunLight;

    public Dropdown marsSatellites;
    public Dropdown jupiterSatellites;

    private RotationAroundPlanet rotationAroundPlaneScript;


    void Start () {
		
        MainCamera = GameObject.FindGameObjectWithTag("MainCamera");
        // marsSatellites.onValueChanged.AddListener(delegate { ValueChangeCheck(marsSatellites); });
        rotationAroundPlaneScript = MainCamera.GetComponent<RotationAroundPlanet>();
    }

    public void ChangeSatellites(string planet)
    {
        switch (planet)
        {
            case "Mars":
                MainCamera.GetComponent<PlanetInfo>().LoadTextToScrollBar(marsSatellites.captionText.text);
                break;
            case "Jupiter":
                MainCamera.GetComponent<PlanetInfo>().LoadTextToScrollBar(jupiterSatellites.captionText.text);
                break;
            default:
                break;
        }
    }
    
   
	public void AssignPlanetCameraCoordinates (string selectedPlanetName)
	{
		planetName = selectedPlanetName;
		
        GameObject planet = GameObject.Find(selectedPlanetName);

        if (planet != null)
        {
            rotationAroundPlaneScript.target = planet.transform;

            

            rotationAroundPlaneScript.distance = baseDistance * planet.transform.localScale.x * planet.transform.parent.localScale.x;
            
           
            rotationAroundPlaneScript.MouseWheelSensitivity = baseScrollSpeed * planet.transform.localScale.x;

            

            rotationAroundPlaneScript.MouseZoomMin = planet.transform.localScale.x * optimalDistance;

            

            Camera.main.fieldOfView = rotationAroundPlaneScript.fovDefault;

            

            if (selectedPlanetName == "Sun")
            {
                SunLight.enabled = false;
            }
            else
            {
                SunLight.enabled = true;
            }
        }
       
	}

}
