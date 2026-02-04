import React, { useState } from 'react';
import { useSession } from '../contexts/SessionContext';
import './SessionManager.css';

interface SessionManagerProps {
  onClose?: () => void;
  onSessionLoadRequested?: (session: { vcfPath?: string; selectedRegion?: string; selectedGenotypes?: string[]; publicVcfFolder?: string; cohortFolder?: string; mode?: string }) => void;
  getSessionData?: () => { vcfPath: string; selectedRegion: string; selectedGenotypes: string[]; publicVcfFolder: string; cohortFolder: string; mode: string };
}

const SessionManager: React.FC<SessionManagerProps> = ({ onClose, onSessionLoadRequested, getSessionData }) => {
  const {
    sessions,
    currentSession,
    saveSession,
    loadSession,
    deleteSession,
    clearCurrentSession,
  } = useSession();

  const [isOpen, setIsOpen] = useState(false);
  const [sessionName, setSessionName] = useState('');
  const [showSaveDialog, setShowSaveDialog] = useState(false);

  const handleSaveCurrent = () => {
    if (!currentSession) {
      setShowSaveDialog(true);
      return;
    }
    saveSession({
      name: currentSession.name,
      vcfPath: currentSession.vcfPath,
      selectedRegion: currentSession.selectedRegion,
      selectedGenotypes: currentSession.selectedGenotypes,
      publicVcfFolder: currentSession.publicVcfFolder,
      cohortFolder: currentSession.cohortFolder,
      mode: currentSession.mode,
      notes: currentSession.notes,
    });
    setShowSaveDialog(false);
    setSessionName('');
  };

  const handleSaveNew = () => {
    if (!sessionName.trim()) {
      alert('Please enter a session name');
      return;
    }
    const data = getSessionData?.() ?? currentSession ?? {
      vcfPath: '',
      selectedRegion: '',
      selectedGenotypes: [] as string[],
      publicVcfFolder: '',
      cohortFolder: '',
      mode: 'individual' as const,
    };
    const mode = (data.mode === 'cohort-read' || data.mode === 'cohort-assembly') ? data.mode : 'individual';
    saveSession({
      name: sessionName,
      vcfPath: data.vcfPath ?? '',
      selectedRegion: data.selectedRegion ?? '',
      selectedGenotypes: data.selectedGenotypes ?? [],
      publicVcfFolder: data.publicVcfFolder ?? '',
      cohortFolder: data.cohortFolder ?? '',
      mode,
    });
    setShowSaveDialog(false);
    setSessionName('');
  };

  const formatDate = (timestamp: number) => {
    return new Date(timestamp).toLocaleString();
  };

  // If onClose is provided, we're being used as a modal (opened from bookmark panel)
  // Otherwise, use the button toggle behavior
  const shouldShowButton = !onClose;
  const isModalOpen = onClose ? true : isOpen;

  return (
    <>
      {shouldShowButton && (
        <button
          className="session-manager-button"
          onClick={() => setIsOpen(!isOpen)}
          title="Manage sessions"
        >
          <span className="session-icon">💾</span>
          <span className="session-label">Sessions</span>
          {currentSession && (
            <span className="session-badge">{currentSession.name}</span>
          )}
        </button>
      )}

      {isModalOpen && (
        <div className="session-manager-overlay" onClick={() => {
          if (onClose) {
            onClose();
          } else {
            setIsOpen(false);
          }
        }}>
          <div className="session-manager-modal" onClick={(e) => e.stopPropagation()}>
            <div className="session-manager-header">
              <h2>💾 Session Manager</h2>
              <button
                className="session-manager-close"
                onClick={() => {
                  if (onClose) {
                    onClose();
                  } else {
                    setIsOpen(false);
                  }
                }}
              >
                ✕
              </button>
            </div>

            <div className="session-manager-content">
              {showSaveDialog ? (
                <div className="session-save-dialog">
                  <h3>Save Session</h3>
                  <input
                    type="text"
                    className="session-name-input"
                    placeholder="Enter session name"
                    value={sessionName}
                    onChange={(e) => setSessionName(e.target.value)}
                    onKeyPress={(e) => {
                      if (e.key === 'Enter') {
                        handleSaveNew();
                      }
                    }}
                    autoFocus
                  />
                  <div className="session-save-actions">
                    <button
                      className="session-btn session-btn-primary"
                      onClick={handleSaveNew}
                    >
                      Save
                    </button>
                    <button
                      className="session-btn session-btn-secondary"
                      onClick={() => {
                        setShowSaveDialog(false);
                        setSessionName('');
                      }}
                    >
                      Cancel
                    </button>
                  </div>
                </div>
              ) : (
                <>
                  <div className="session-manager-actions">
                    <button
                      className="session-btn session-btn-primary"
                      onClick={handleSaveCurrent}
                    >
                      💾 Save Current Session
                    </button>
                    {currentSession && (
                      <button
                        className="session-btn session-btn-secondary"
                        onClick={clearCurrentSession}
                      >
                        Clear Current
                      </button>
                    )}
                  </div>

                  {currentSession && (
                    <div className="current-session-info">
                      <h3>Current Session</h3>
                      <div className="session-info-item">
                        <strong>Name:</strong> {currentSession.name}
                      </div>
                      <div className="session-info-item">
                        <strong>Saved:</strong> {formatDate(currentSession.timestamp)}
                      </div>
                      {currentSession.vcfPath && (
                        <div className="session-info-item">
                          <strong>VCF:</strong> {currentSession.vcfPath}
                        </div>
                      )}
                      {currentSession.selectedRegion && (
                        <div className="session-info-item">
                          <strong>Region:</strong> {currentSession.selectedRegion}
                        </div>
                      )}
                    </div>
                  )}

                  <div className="session-list">
                    <h3>Saved Sessions ({sessions.length})</h3>
                    {sessions.length === 0 ? (
                      <p className="session-empty">No saved sessions</p>
                    ) : (
                      <div className="session-items">
                        {sessions.map((session) => (
                          <div
                            key={session.id}
                            className={`session-item ${
                              currentSession?.id === session.id ? 'active' : ''
                            }`}
                          >
                            <div className="session-item-content">
                              <div className="session-item-name">{session.name}</div>
                              <div className="session-item-meta">
                                {formatDate(session.timestamp)}
                              </div>
                              {session.selectedRegion && (
                                <div className="session-item-region">
                                  Region: {session.selectedRegion}
                                </div>
                              )}
                            </div>
                            <div className="session-item-actions">
                              <button
                                className="session-item-btn"
                                onClick={() => {
                                  loadSession(session.id);
                                  onSessionLoadRequested?.(session);
                                  onClose?.();
                                }}
                                title="Load session"
                              >
                                📂
                              </button>
                              <button
                                className="session-item-btn session-item-btn-danger"
                                onClick={() => {
                                  if (window.confirm(`Delete session "${session.name}"?`)) {
                                    deleteSession(session.id);
                                  }
                                }}
                                title="Delete session"
                              >
                                🗑️
                              </button>
                            </div>
                          </div>
                        ))}
                      </div>
                    )}
                  </div>
                </>
              )}
            </div>
          </div>
        </div>
      )}
    </>
  );
};

export default SessionManager;




