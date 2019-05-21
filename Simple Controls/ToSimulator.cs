using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class ToSimulator : MonoBehaviour {

    public void ChangeScene(string sceneName)
    {
        //chama a cena nova
        Application.LoadLevel (sceneName);
    }
}
