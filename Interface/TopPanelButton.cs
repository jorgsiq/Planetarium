using UnityEngine;

public class TopPanelButton : MonoBehaviour {

    private Animator buttonAnimator;

    void Start()
    {
        buttonAnimator = this.GetComponent<Animator>();
    }

    public void HoverButton()
    {
        if (buttonAnimator.GetCurrentAnimatorStateInfo(0).IsName("TB Hover to Pressed"))
        {
 
        }

        else
        {
            buttonAnimator.Play("TB Hover");
        }
    }

    public void NormalizeButton()
    {
        if (buttonAnimator.GetCurrentAnimatorStateInfo(0).IsName("TB Hover to Pressed"))
        {

        }

        else
        {
            buttonAnimator.Play("TB Hover to Normal");
        }
    }
}
