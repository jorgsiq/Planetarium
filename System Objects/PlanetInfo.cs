using UnityEngine;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Text;
using UnityEngine.UI;

public class PlanetInfo : MonoBehaviour {

    private List<string> infoList { get; set; }
    public Text contentText;
    public Text planetNameText;
    private RectTransform rectTransform;
    public float dpi = 5;
    private PlanetSwitch planetSwitchScript;



	void Start () {
        rectTransform = contentText.GetComponent<RectTransform>();
        planetSwitchScript = gameObject.AddComponent<PlanetSwitch>();
        planetSwitchScript.SunLight = GameObject.Find("SunLight").GetComponent<Light>();
    }

    public void LoadTextToScrollBar(string name)
    {
        object resourceFile = Resources.Load(name);
        
    
        contentText.text = string.Empty;

        if(resourceFile != null)
        {
           
            resourceFile = resourceFile.ToString().Replace("TAB", "\t");
            contentText.text += resourceFile;
            
            rectTransform.sizeDelta = new Vector2(rectTransform.sizeDelta.x, resourceFile.ToString().Length * dpi);
        }
        else
        {
         
            contentText.text = string.Format("Please add {0}.txt to Resources folder", name);
        }

        planetNameText.text = name;
        planetSwitchScript.AssignPlanetCameraCoordinates(name);
    }
}
